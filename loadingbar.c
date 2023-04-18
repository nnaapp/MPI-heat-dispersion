#include <stdio.h>
#include <stdlib.h>
#include "loadingbar.h"

// initialized loading bars (maximum length of 512)
// Basically a constructor, but this is C and I wouldn't make it a class anyway.
struct LoadingBar loadingbar_init(int max, char Fill, char Empty, char left, char right)
{
    struct LoadingBar bar;

    bar.curLen = 0;
    bar.maxLen = max;
    if (bar.maxLen > 512)
        bar.maxLen = 512;
    bar.percent = 0;

    bar.fill = Fill;
    bar.empty = Empty;

    bar.leftCap = left;
    bar.rightCap = right;

    return bar;
}

// draws loading bars
// Draws on the same line, using \r to return to start,
// so as to not be ugly and clutter the command prompt.
// Buffer must be flushed, or it will not work right due to
// the usage of \r, I believe.
void loadingbar_draw(struct LoadingBar *bar)
{
    // move to start of line and print left cap char
    if (bar->curLen > bar->maxLen)
        bar->curLen = bar->maxLen;

    if (bar->percent > 100)
        bar->percent = 100;

    char *buffer = calloc(bar->maxLen, 1);
    int pos = 0;
    for (int i = 0; i < bar->curLen; i++)
    {
        buffer[i] = bar->fill;
        pos++;
    }

    for (int i = pos; i < bar->maxLen; i++)
    {
        buffer[i] = bar->empty;
        pos++;
    }

    // print breakdown:
    // \33[2K - clear current line, works on MOST consoles (if this breaks, its probably your console)
    // \r     - carriage return, goes back to start of current line
    // %c%d%% - left cap character, then percentage completion, then % sign escaped to a % literal
    // %s%c   - string buffer, avoids printing 50 different characters one by one, and then right cap char
    printf("\33[2K\r%c%d%% %s%c", bar->leftCap, bar->percent, buffer, bar->rightCap);

    fflush(stdout);
}