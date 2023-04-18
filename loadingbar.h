#ifndef LOADING_BAR_H
#define LOADING_BAR_H

struct LoadingBar
{
    int maxLen, curLen, percent;
    char fill, empty;
    char leftCap, rightCap;
};

struct LoadingBar loadingbar_init(int, char, char, char, char);
void loadingbar_draw(struct LoadingBar *);

#endif