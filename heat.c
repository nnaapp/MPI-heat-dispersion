#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "loadingbar.h" // defines and draws progress bar in console, gives user something to stare at
#include "heat.h"

void handle_loading_bar(int, int, struct LoadingBar *);

void fill_heaters(float *, struct Heater *, int, int, int, int, int, int);
void get_heater_count(int *, char *, int, int);
void get_heaters(struct Heater **, char *, int, int, int);

void matrix_distribute(float *, float *, int, int, int, int, int, int);
void matrix_stitch(float *, float *, int, int, int, int, int, int);
void matrix_halo(float *, float *, float *, int, int, int, int, int, int, float);

float *matrix_init_empty(int, int);
float *matrix_init(int, int, float);
void matrix_out(float *, int, int, char *);
void matrix_step(float **, float**, float *, float *, int, int, float, float, int, int);
float matrix_sum_neighbors(float *, float *, float *, int, int, int, int, float);

void mpi_try_exit();

int main(int argc, char **argv)
{
    if (argc != EXPECTED_ARGS)
    {
        printf("Invalid usage.\n");
        printf("Example: ./heat numRows numCols k basetemp timesteps heaterFileName outputFileName\n");

        mpi_try_exit();
        return 1;
    }

    /* Command line arguments and parsing*/
    int numRows, numCols, timesteps;
    float baseTemp, transferRate;
    char *heaterFileName;
    char *outFileName;
    char *outImgName;

    // these are done on every process
    numRows = atoi(argv[1]);
    numCols = atoi(argv[2]);
    timesteps = atoi(argv[5]);
    char *ptr; // ptr for strtod
    baseTemp = strtod(argv[4], &ptr);
    transferRate = strtod(argv[3], &ptr);

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int normProcHeight = numRows / size;
    int remainder = numRows - (normProcHeight * size);
    int localHeight = normProcHeight;

    if (rank == size - 1)
        localHeight += remainder;

    // these are done on the main process,
    // and are done after init because we need to know what process we are on
    if (rank == 0)
    {
        heaterFileName = argv[6];
        outFileName = argv[7];
        outImgName = (char *)malloc(sizeof(char) * strlen(outFileName) + 4);
        strcpy(outImgName, outFileName);
        strcat(outImgName, ".bmp");
    }


    /* Argument validation and error prevention */
    if (numRows < 1 || numCols < 1)
    {
        printf("Invalid matrix dimensions, must be 1x1 or greater.\n");
        mpi_try_exit();
        return 1;
    }

    if (timesteps < 1)
    {
        printf("Invalid number of timesteps, must be >0, time can't go backwards.\n");
        mpi_try_exit();
        return 1;
    }

    if (transferRate > TRANSFER_MAX || transferRate < TRASNFER_MIN)
    {
        printf("Invalid heat transfer rate, choose a number between 1 and 1.1 (inclusive).\n");
        mpi_try_exit();
        return 1;
    }


    /* File reading, data and variable initialization */

    // heaters is a pointer array of heater structs to contain the row/col/temp of each heater
    // dynamic arrays are impossible to calculate the length of in code, so that is stored from file
    int heaterCount = 0;
    struct Heater *heaters = NULL;
    get_heater_count(&heaterCount, heaterFileName, rank, size);
    get_heaters(&heaters, heaterFileName, heaterCount, rank, size);

    // initialize matrix of argument size and temp, fill it with heaters from file
    // these matrices persist, and are swapped around instead of re-allocated
    float *matrix = NULL;
    struct LoadingBar progress;
    if (rank == 0)
    {
        matrix = matrix_init(numCols, numRows, baseTemp);
        fill_heaters(matrix, heaters, heaterCount, numCols, numRows, numRows, 0, 0);

        progress = loadingbar_init(50, '#', '-', '[', ']');
        loadingbar_draw(&progress);
    }

    float *localMatrix = malloc(sizeof(*localMatrix) * (localHeight * numCols));
    float *gTop = malloc(sizeof(*gTop) * numCols);
    float *gBottom = malloc(sizeof(*gBottom) * numCols);
    float *localTemp = malloc(sizeof(*localTemp) * (localHeight * numCols));
    
    matrix_distribute(localMatrix, matrix, rank, size, normProcHeight, remainder, localHeight, numCols);
    for (int i = 0; i < timesteps; i++)
    {
        matrix_halo(localMatrix, gBottom, gTop, rank, size, normProcHeight, remainder, localHeight, numCols, baseTemp);
        matrix_step(&localMatrix, &localTemp, gTop, gBottom, numCols, localHeight, transferRate, baseTemp, rank, size);

        fill_heaters(localMatrix, heaters, heaterCount, numCols, normProcHeight, localHeight, rank, size);

        if (rank == 0)
        {
            handle_loading_bar(i, timesteps - 1, &progress);
        }
    }
    matrix_stitch(localMatrix, matrix, rank, size, normProcHeight, remainder, localHeight, numCols);

    if (rank == 0)
    {
        printf("\n");
        matrix_out(matrix, numCols, numRows, outFileName);
        free(matrix);
        free(heaters);
    }

    free(localMatrix);
    free(localTemp);
    free(gTop);
    free(gBottom);

    /* Matrix timesteps, data processing into CSV and BMP image */
    MPI_Finalize();
    return 0;
}


void mpi_try_exit()
{
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized)
        MPI_Finalize();
}

// Takes a 2d array matrix, array of heaters, and the number of heaters.
// Returns the matrix with the heaters placed where they belong, based on struct.
void fill_heaters(float *matrix, struct Heater *heaters, int arrayLen, int cols, int localHeight, int normHeight, int rank, int size)
{
    int start = (rank * normHeight) * cols;
    int end = (start + (localHeight * cols)) - 1;
    int span = end - start;
    for (int i = 0; i < arrayLen; i++)
    {
        int relative_pos = (heaters[i].col + (heaters[i].row * cols)) - start;
        if (relative_pos > 0 && relative_pos < span)
            matrix[relative_pos] = heaters[i].temp;
    }
}

// Reads the first line of the named file into a buffer,
// the first line is always the number of heaters,
// and returns that number as an int.
void get_heater_count(int *heaterCount, char *heaterFileName, int rank, int size)
{
    if (rank == 0)
    {
        FILE *heaterFile;
        heaterFile = fopen(heaterFileName, "r");

        if (!heaterFile)
        {
            *heaterCount = 0;
            for (int i = 0; i < size; i++)
            {
                MPI_Send(heaterCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            return;
        }

        char buffer[READ_BUFFER]; // byte buffer to store the first line in

        fgets(buffer, READ_BUFFER, heaterFile); // first line (by new line char) contains number of heaters
        *heaterCount = atoi(buffer);   // aka number of lines to read
        fclose(heaterFile);

        for (int i = 0; i < size; i++)
        {
            MPI_Send(heaterCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        return;
    }

    MPI_Recv(heaterCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// Takes file name and number of heaters,
// returns null if no heaters, and parses each line of the file
// given by the number of heaters into ints coordinates and float temp.
// Returns heater struct with heater data.
void get_heaters(struct Heater **heaters, char *heaterFileName, int numHeaters, int rank, int size)
{
    if (numHeaters == 0)
    {
        *heaters = NULL;
        return;
    }

    *heaters = malloc(sizeof(struct Heater) * numHeaters);

    if (rank == 0)
    {
        FILE *heaterFile;
        heaterFile = fopen(heaterFileName, "r"); // r for read only
        char buffer[READ_BUFFER];

        fscanf(heaterFile, "%s", buffer);
        for (int i = 0; i < numHeaters; i++)
        {
            fscanf(heaterFile, "%s", buffer);
            int row = atoi(buffer);

            fscanf(heaterFile, "%s", buffer);
            int col = atoi(buffer);

            fscanf(heaterFile, "%s", buffer);
            char *ptr;
            float t = strtod(buffer, &ptr);

            for (int i = 0; i < size; i++)
            {
                MPI_Send(&row, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&col, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&t, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }

            struct Heater tmp;
            tmp.row = row;
            tmp.col = col;
            tmp.temp = t;

            (*heaters)[i] = tmp;
        }

        fclose(heaterFile);
        return;
    }

    for (int i = 0; i < numHeaters; i++)
    {
        int row = 0, col = 0;
        float t = 0;

        MPI_Recv(&row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&col, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        struct Heater tmp;
        tmp.row = row;
        tmp.col = col;
        tmp.temp = t;

        (*heaters)[i] = tmp;
    }
}

void matrix_distribute(float *localMatrix, float *matrix, int rank, int size, int normHeight, int heightRem, int localHeight, int cols)
{
    if (rank == 0)
    {
        for (int i = 0; i < normHeight * cols; i++)
        {
            localMatrix[i] = matrix[i];
        }

        if (size - 1 > 0)
        {
            for (int i = 1; i < size; i++)
            {
                int curHeight = (i == size - 1) ? normHeight + heightRem : normHeight;
                int area = curHeight * cols;
                float *tmp = matrix + ((normHeight * i) * cols);
                MPI_Send(tmp, area, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
        }
    }

    if (rank != 0)
    {
        int area = localHeight * cols;
        MPI_Recv(localMatrix, area, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void matrix_stitch(float *localMatrix, float *matrix, int rank, int size, int normHeight, int heightRem, int localHeight, int cols)
{
    if (rank == 0)
    {
        int localArea = localHeight * cols;
        for (int i = 0; i < localArea; i++)
        {
            matrix[i] = localMatrix[i];
        }

        for (int i = 1; i < size; i++)
        {
            int curHeight = (i == size - 1) ? normHeight + heightRem : normHeight;
            int area = curHeight * cols;
            float *tmp = matrix + ((normHeight * i) * cols);
            MPI_Recv(tmp, area, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (rank != 0)
    {
        int area = localHeight * cols;
        MPI_Send(localMatrix, area, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }
}

void matrix_halo(float *localMatrix, float *gBottom, float *gTop, int rank, int size, int normHeight, int heightRem, int localHeight, int cols, float base)
{
    MPI_Request sendTop, sendBottom, recTop, recBottom;
    if (rank == 0)
    {
        for (int i = 0; i < cols; i++)
        {
            gTop[i] = base;
        }

        float *tmp = localMatrix + ((localHeight - 1) * cols);

        if (size > 1)
        {
            MPI_Isend(tmp, cols, MPI_FLOAT, (rank + 1) % size, 0, MPI_COMM_WORLD, &sendBottom);
            MPI_Irecv(gBottom, cols, MPI_FLOAT, (rank + 1) % size, 0, MPI_COMM_WORLD, &recBottom);
            MPI_Wait(&sendBottom, MPI_STATUS_IGNORE);
            MPI_Wait(&recBottom, MPI_STATUS_IGNORE);
        }
    }

    if (rank == size - 1)
    {
        for (int i = 0; i < cols; i++)
        {
            gBottom[i] = base;
        }

        float *tmp = localMatrix;

        if (size > 1)
        {
            MPI_Isend(tmp, cols, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &sendTop);
            MPI_Irecv(gTop, cols, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &recTop);
            MPI_Wait(&sendTop, MPI_STATUS_IGNORE);
            MPI_Wait(&recTop, MPI_STATUS_IGNORE);
        }
    }

    if (rank > 0 && rank < size - 1)
    {
        float *tmpTop = localMatrix;
        float *tmpBottom = localMatrix + ((localHeight - 1) * cols);

        MPI_Isend(tmpTop, cols, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &sendTop);
        MPI_Irecv(gTop, cols, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &recTop);
        MPI_Wait(&sendTop, MPI_STATUS_IGNORE);
        MPI_Wait(&recTop, MPI_STATUS_IGNORE);
        MPI_Isend(tmpBottom, cols, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &sendBottom);
        MPI_Irecv(gBottom, cols, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &recBottom);
        MPI_Wait(&sendBottom, MPI_STATUS_IGNORE);
        MPI_Wait(&recBottom, MPI_STATUS_IGNORE);
    }
}

// Takes row/col sizes, allocates an EMPTY matrix accordingly.
// Returns matrix ptr
float *matrix_init_empty(int cols, int rows)
{
    // 1d array which will be indexed like a 2d array
    float *matrix_ptr = (float *)malloc((cols * rows) * sizeof(float));

    return matrix_ptr;
}

// Takes row/col sizes and the base temp, and allocates a matrix accordingly.
// Returns matrix ptr
float *matrix_init(int cols, int rows, float base)
{
    // each column is a pointer to an array (pointer)
    float *matrix_ptr = (float *)malloc((cols * rows) * sizeof(float));

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            matrix_ptr[j + (i * cols)] = base; // base is default temperature
        }
    }

    return matrix_ptr;
}

// Takes 2d array matrix, its dimensions, and an output file name.
// Prints a buffer containing every cell of the matrix to the given file.
// Done this way to avoid literally 25 million fprintf's, because thats slow.
void matrix_out(float *matrix, int cols, int rows, char *outFileName)
{
    int writeBuffSize = (cols * rows) * (sizeof(char) * WRITE_BUFF_MULT);
    char *write_buffer = (char *)malloc(writeBuffSize);
    write_buffer[0] = '\0';
    char convert_buffer[CONV_BUFF_SIZE];
    int strSize = 0;

    FILE *outFile;
    outFile = fopen(outFileName, "w"); // w for write only

    if (!outFile) // !outFile = file could not be opened, is NULL
    {
        printf("ERROR: Output file could not be opened.\n");
        return;
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            // Converts float to string and stores it in a 64 byte buffer
            snprintf(convert_buffer, sizeof(convert_buffer), "%.1f,", matrix[j + (i * cols)]);
            // Concatonates to main buffer by null terminator
            int numLen = strlen(convert_buffer);
            for (int k = 0; k < numLen; k++)
            {
                write_buffer[k + strSize] = convert_buffer[k];
            }
            write_buffer[numLen + strSize] = '\0';
            strSize += numLen;
        }

        // Checks if buffer is 80% full, and writes/empties if so
        if (strSize > (writeBuffSize * 0.8))
        {
            fprintf(outFile, "%s", write_buffer);
            memset(write_buffer, 0, strSize);
            write_buffer[0] = '\0';
            strSize = 0;
        }
        write_buffer[strSize] = '\n';
        strSize++;
    }
    // ONE write to file, because we are going for speed here
    fprintf(outFile, "%s", write_buffer);

    free(write_buffer);
}

// Takes ADDRESS of matrix (this is necessary for efficient swapping and avoiding memory leaks)
// as well as dimensions of matrix, transfer rate, temperature, and thread count.
// Performs one time step on the array using given temp/rate/dimensions.
void matrix_step(float **matrix, float **tmpMatrix, float *gTop, float *gBottom, int cols, int localHeight, float k, float base, int rank, int size)
{
    float *newMatrix = *tmpMatrix;
    float *curMatrix = *matrix; // derefence address of matrix to usable form

    // Each thread is given one "cell" (x/y coordinate) at a time
    // and calculates the new temperature based on neighbors.
    // This new value is stored in a temporary matrix, and doing this
    // requires no writes to the original matrix. This means
    // there are no race conditions here, as no thread will write
    // where any other threat wants to write.
    for (int i = 0; i < localHeight; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            newMatrix[j + (i * cols)] = 
                (curMatrix[j + (i * cols)] + (k * matrix_sum_neighbors(curMatrix, gTop, gBottom, j, i, cols, localHeight, base)) / 8.0) / 2.0;
        }
    }

    float *tmp = *matrix;
    *matrix = *tmpMatrix; // put tmpMatrix at the address of main matrix
    *tmpMatrix = tmp;
}

// Takes a matrix (float *), coordinates, dimensions, and a default temperature.
// Does not need parallelized as the smallest reasonable chunks each thread can do
// are to calculate each index's neighbor sum.
// Returns the sum of the coordinate's neighbors, defaulting out-of-bounds indices
// to the default temperature.
float matrix_sum_neighbors(float *matrix, float *gTop, float *gBottom, int x, int y, int cols, int localHeight, float base)
{
    float sum = 0;

    // sort of ugly, but makes the common case quite fast.
    // very little comparisons/math overhead, just raw summing.
    if (x > 0 && x < cols - 1 && y > 0 && y < localHeight - 1)
    {
        sum += matrix[x - 1 + ((y - 1) * cols)];
        sum += matrix[x     + ((y - 1) * cols)];
        sum += matrix[x - 1 + ( y      * cols)];

        sum += matrix[x + 1 + ((y + 1) * cols)];
        sum += matrix[x     + ((y + 1) * cols)];
        sum += matrix[x + 1 + ( y      * cols)];

        sum += matrix[x - 1 + ((y + 1) * cols)];
        sum += matrix[x + 1 + ((y - 1) * cols)];

        return sum;
    }

    // if above comparison failed, we are on some sort of edge.
    // this loop will only be encountered a number of times equal to the
    // perimeter of the matrix.
    // basically just goes over the neighbors and picks out the bad ones.
    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            if (i == 0 && j == 0)
                continue;

            int cur_x = x + j;
            int cur_y = y + i;

            // If the current Y is out of bounds, but not the current X
            if ((cur_y < 0 || cur_y >= localHeight) && (cur_x >= 0 && cur_x < cols))
            {
                sum += (cur_y < 0) ? gTop[cur_x] : gBottom[cur_x];
                continue;
            }

            // If the current coordinates are otherwise out of bounds, just add the base
            if (cur_x < 0 || cur_x >= cols || cur_y < 0 || cur_y >= localHeight)
            {
                sum += base;
                continue;
            }

            sum += matrix[cur_x + (cur_y * cols)];
        }
    }

    return sum;
}

// Handles loading bar, checks if it needs an update.
// Conditions for update are an increase in whole-number percent,
// or another filling-character needing to be placed.
void handle_loading_bar(int loop_pos, int timesteps, struct LoadingBar *bar)
{
    unsigned long curProgress = (((unsigned long)loop_pos) * bar->maxLen) / timesteps;
    int curPercent = floor(((float)(loop_pos)) / ((float)(timesteps)) * 100);

    if (curProgress != bar->curLen || curPercent > bar->percent)
    {
        bar->curLen = curProgress;
        bar->percent = curPercent;
        loadingbar_draw(bar);
    }
}
