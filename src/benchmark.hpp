#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <ctime>
#include <cmath>

double static readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}




#endif