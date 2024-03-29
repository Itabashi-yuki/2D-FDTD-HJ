#include "fdtd2d.h"
#include <fstream>
#include <iostream>

void judge_Ez_div(double **Ez, int n, int n0, double *Ez_max, double Ne, double nu){
    if(n % 5 == 0){
        for(int i = Nx_iono_lower; i <= Nx_iono_upper; i+=2){
            if(*Ez_max < abs(Ez[i][Ny/2])){
                *Ez_max = abs(Ez[i][Ny/2]);
            }
        }
        std::ofstream ofs_max("./data/" + global_dirName + "/Ez_max_dt="  + "_Ne=" + std::to_string(Ne) + "_Tmax=" + std::to_string(Tmax) + ".dat",std::ios::app);
        ofs_max << dt * n << " " << *Ez_max << std::endl;
        ofs_max.close();
    }
}