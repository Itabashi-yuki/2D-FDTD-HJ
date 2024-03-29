#include "fdtd2d.h"
#include <iostream>
#include <fstream>

void output_obs(double **Ex, double **Ey, double **Ez, double **Hx, double **Hy, double **Hz,
            double ***Jex, double ***Jey, double ***Jez, int n){
    std::ofstream ofs("./data/" + global_dirName + "/obs.dat",std::ios::app);
    std::ofstream ofs_2("./data/" + global_dirName + "/obs2.dat",std::ios::app);

    ofs << n * dt << " " << Ex[int(obs_x / dx)][int(obs_y / dy)] << " " << Ey[int(obs_x / dx)][int(obs_y / dy)] << " " <<
         Ez[int(obs_x / dx)][int(obs_y / dy)] << " " << Hx[int(obs_x / dx)][int(obs_y / dy)] << " " <<
          Hy[int(obs_x / dx)][int(obs_y / dy)] << " " << Hz[int(obs_x / dx)][int(obs_y / dy)] << " " << std::endl;

    ofs_2 << n * dt << " " << Ex[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
            << Ey[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
            << Ez[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
            << Hx[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
            << Hy[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
            << Hz[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " << std::endl;

    ofs.close();
    ofs_2.close();
}