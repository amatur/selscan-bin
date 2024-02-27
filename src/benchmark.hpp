#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <ctime>
#include <cmath>

double readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


#endif