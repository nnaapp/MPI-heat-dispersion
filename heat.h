#ifndef HEAT_H
#define HEAT_H

#define EXPECTED_ARGS 8
#define TRANSFER_MAX 1.1000001 // floating point imprecision, man
#define TRASNFER_MIN 1

#define READ_BUFFER 1024

#define WRITE_BUFF_MULT 8
#define CONV_BUFF_SIZE 64

#define IMG_DEFAULT 1024
#define IMG_MAX IMG_DEFAULT*5
#define BPP 3
                                // best results:
#define B_DEFAULT 10            // 10
#define G_DEFAULT 127           // 127
#define G_DEFAULT_BRIGHT 128    // 128
#define R_DEFAULT 10            // 10
                                //
#define RED_SHIFT 3             // 3
#define BLUE_SHIFT -5           // -4
#define GREEN_SHIFT 2           // 2

struct Heater
{
    int row;
    int col;
    float temp;
};

#endif