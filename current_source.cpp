#include <cmath>
#include "fdtd2d.h"

double Jz(double t){
    return std::exp( - (t - t0) * (t - t0) / 2.0 / sigma / sigma );

    // double Jz = 0.0;
    // if(0 <= t && t < t0){
    //     Jz = sin( 2 * M_PI * f0 * t) * std::exp( -( t - t0) * ( t - t0) / 2.0 / sigma / sigma) / std::sqrt( dx * dx + dy * dy);
    // } else if (t0 < t){
    //     Jz = sin(2 * M_PI * f0 * t) / std::sqrt( dx * dx + dy * dy);
    // }
    // return Jz;

}