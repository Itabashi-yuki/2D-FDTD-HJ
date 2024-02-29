#include <cmath>
#include "fdtd2d.h"

void current_source(double *Jz, int n){
    double t = n * dt;
    double sigma = 15 * dt;
    double t0 = 6 * sigma;

    Jz[n] = std::exp( - (t - t0) * (t - t0) / 2.0 / sigma / sigma );

}