#include "fdtd2d.h"
#include <iostream>
#include <fstream>

void output_obs(double **Ex, double **Ey, double **Ez, double **Hx, double **Hy, double **Hz,
            double ***Jex, double ***Jey, double ***Jez, int n){
    std::ofstream ofs("./data/" + global_dirName + "/obs.dat",std::ios::app);
    std::ofstream ofs_2("./data/" + global_dirName + "/obs2.dat",std::ios::app);
    std::ofstream ofs_3("./data/" + global_dirName + "/obs3.dat",std::ios::app);
    std::ofstream ofs_4("./data/" + global_dirName + "/obs4.dat",std::ios::app);

    ofs << "# (70,50)における時間変化" << std::endl;
    ofs << n * dt << " " << Ex[int(obs_x / dx)][int(obs_y / dy)] << " " << Ey[int(obs_x / dx)][int(obs_y / dy)] << " " <<
         Ez[int(obs_x / dx)][int(obs_y / dy)] << " " << Hx[int(obs_x / dx)][int(obs_y / dy)] << " " <<
          Hy[int(obs_x / dx)][int(obs_y / dy)] << " " << Hz[int(obs_x / dx)][int(obs_y / dy)] << " " << std::endl;

//     ofs_2 << n * dt << " " << Ex[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
//             << Ey[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
//             << Ez[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
//             << Hx[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
//             << Hy[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " 
//             << Hz[int(source_x / dx) - 10][int(source_y / dy) - 10] << " " << std::endl;
    ofs_2 << "# (70,70)における時間変化" << std::endl;
    ofs_2 << n * dt << " " << Ex[int(obs_x / dx)][int(obs_y / dy) + 40] << " " << Ey[int(obs_x / dx)][int(obs_y / dy) + 40] << " " <<
         Ez[int(obs_x / dx)][int(obs_y / dy) + 40] << " " << Hx[int(obs_x / dx)][int(obs_y / dy) + 40] << " " <<
          Hy[int(obs_x / dx)][int(obs_y / dy) + 40] << " " << Hz[int(obs_x / dx)][int(obs_y / dy) + 40] << " " << std::endl;

    ofs_3 << "# (70,30)における時間変化" << std::endl;
    ofs_3 << n * dt << " " << Ex[int(obs_x / dx)][int(obs_y / dy) - 40] << " " << Ey[int(obs_x / dx)][int(obs_y / dy) - 40] << " " <<
         Ez[int(obs_x / dx)][int(obs_y / dy) - 40] << " " << Hx[int(obs_x / dx)][int(obs_y / dy) - 40] << " " <<
          Hy[int(obs_x / dx)][int(obs_y / dy) - 40] << " " << Hz[int(obs_x / dx)][int(obs_y / dy) - 40] << " " << std::endl;

    ofs_4 << "# (50,90)における時間変化" << std::endl;
    ofs_4 << n * dt << " " << Ex[int(obs_x / dx) + 40][int(obs_y / dy)] << " " << Ey[int(obs_x / dx) + 40][int(obs_y / dy)] << " " <<
         Ez[int(obs_x / dx) + 40][int(obs_y / dy)] << " " << Hx[int(obs_x / dx) + 40][int(obs_y / dy)] << " " <<
          Hy[int(obs_x / dx) + 40][int(obs_y / dy)] << " " << Hz[int(obs_x / dx) + 40][int(obs_y / dy)] << " " << std::endl;

    ofs.close();
    ofs_2.close();
    ofs_3.close();
    ofs_4.close();
}