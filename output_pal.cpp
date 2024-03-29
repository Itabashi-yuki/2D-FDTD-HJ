#include "fdtd2d.h"
#include <iostream>
#include <fstream>
#include <string>

void output_pal(double Ne, double nu){
    std::ofstream ofs_pal("./data/" + global_dirName + "/pal.dat",std::ios::app);

    ofs_pal << "(pi / 2, 0)から(pi / 4, pi, 4)に向きを変更" << std::endl;
    ofs_pal << std::endl;

    ofs_pal << "Rx = " << Rx << std::endl;
    ofs_pal << "Ry = " << Ry << std::endl;
    ofs_pal << "dx = " << dx << std::endl;
    ofs_pal << "dy = " << dy << std::endl;
    ofs_pal << "sorce_x = " << source_x << std::endl;
    ofs_pal << "sorce_y = " << source_y << std::endl;
    ofs_pal << "obs_x = " << obs_x << std::endl;
    ofs_pal << "obs_y = " << obs_y << std::endl;
    ofs_pal << "dt = " << dt << std::endl;
    ofs_pal << "Tmax = " << Tmax << std::endl;
    ofs_pal << "xi = " << xi << std::endl;
    // ofs_pal << "Cdt = " << Cdt << std::endl;

    ofs_pal << std::endl;
    ofs_pal << "---PML パラメタ---" << std::endl;
    ofs_pal << "M = " << M << std::endl;
    ofs_pal << "R = " << R << std::endl;
    ofs_pal << "L = " << L << std::endl;

    ofs_pal << std::endl;
    ofs_pal << "---Plasma パラメタ---" << std::endl;
    ofs_pal << "Rx_iono_lower = " << Rx_iono_lower << std::endl;
    ofs_pal << "Rx_iono_upper = " << Rx_iono_upper << std::endl;
    ofs_pal << "Ry_iono_lower = " << Ry_iono_lower << std::endl;
    ofs_pal << "Ry_iono_upper = " << Ry_iono_upper << std::endl;
    ofs_pal << "B0 = " << B0 << std::endl;
    ofs_pal << "THETA = " << std::to_string(THETA) << std::endl;
    ofs_pal << "PHI = " << std::to_string(PHI) << std::endl;
    ofs_pal << "Ne = " << cal_Ne(Ne) << std::endl;
    ofs_pal << "nu = " << cal_nu(nu) << std::endl;

    ofs_pal.close();
}