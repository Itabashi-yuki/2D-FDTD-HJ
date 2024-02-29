#include "fdtd2d.h"
#include <fstream>
#include <iostream>

void initialize_PML(double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2, double *CHX1, double *CHX2, double *CHY1, double *CHY2){
    double sigma_max = - (M + 1.0) * Y / 2.0 / L / dx * std::log(R);
    double *sigma_Ezx = new double[Nx];
    double *sigma_Ezy = new double[Nx];
    double *sigma_Hx = new double[Ny];
    double *sigma_Hy = new double[Nx];

    /* σの設定 */
    for(int i = 1; i <= Nx - 1; i++){
        sigma_Ezx[i] = 0.0;
    }

    for(int j = 1; j <= Ny - 1; j++){
        sigma_Ezy[j] = 0.0;
    }

    for(int j = 0; j < Ny; j++){
        sigma_Hx[j] = 0.0;
    }

    for(int i = 0; i < Nx; i++){
        sigma_Hy[i] = 0.0;
    }

    for(int i = 1; i <= L; i++){
        sigma_Ezx[i] = sigma_max * std::pow(((L*dx - i*dx) / L / dx), M);
    }

    for(int j = 1; j <= L; j++){
        sigma_Ezy[j] = sigma_max * std::pow(((L*dy - j*dy) / L / dy), M);
    }

    for(int j = 0; j < L; j++){
        sigma_Hx[j] = Z * Z * sigma_max * std::pow(((L*dy - (j + 0.5) * dy) / L / dy), M);
    }

    for(int i = 0; i < L; i++){
        sigma_Hy[i] = Z * Z * sigma_max * std::pow(((L*dx - (i + 0.5)*dx) / L / dx), M);
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        sigma_Ezx[i] = sigma_max * std::pow(((i*dx - (Nx - L) * dx) / L / dx), M);
    }

    for(int j = Ny - L; j <= Ny - 1; j++){
        sigma_Ezy[j] = sigma_max * std::pow(((j*dy - (Ny - L) * dy) / L / dy), M);
    }

    for(int j = Ny - L; j < Ny; j++){
        sigma_Hx[j] = Z * Z * sigma_max * std::pow((((j + 0.5) * dy - (Ny - L) * dy) / L / dy), M);
    }

    for(int i = Nx - L; i < Nx; i++){
        sigma_Hy[i] = Z * Z * sigma_max * std::pow((((i + 0.5) * dx - (Nx - L) * dx) / L / dx), M);
    }


    std::ofstream ofs_CEZX("CEZX.dat");
    std::ofstream ofs_CEZY("CEZY.dat");
    std::ofstream ofs_CHX("CHX.dat");
    std::ofstream ofs_CHY("CHY.dat");

    for(int i = 1; i <= Nx - 1; i++){
        CEZX1[i] = (EPS0 / dt - sigma_Ezx[i] / 2.0) / (EPS0 / dt + sigma_Ezx[i] / 2.0);
        CEZX2[i] = 1.0 / (EPS0 / dt + sigma_Ezx[i] / 2.0) / dx;

        ofs_CEZX << i << " " << CEZX1[i] << " " << CEZX2[i] << " " << sigma_Ezx[i] << std::endl;
    }

    for(int j = 1; j <= Ny - 1; j++){
        CEZY1[j] = (EPS0 / dt - sigma_Ezy[j] / 2.0) / (EPS0 / dt + sigma_Ezy[j] / 2.0);
        CEZY2[j] = 1.0 / (EPS0 / dt + sigma_Ezy[j] / 2.0) / dy;

        ofs_CEZY << j << " " << CEZX1[j] << " " << CEZX2[j] << " " << sigma_Ezy[j] << std::endl;
    }

    for(int j = 0; j < Ny; j++){
        CHX1[j] = (MU0 / dt - sigma_Hx[j] / 2.0) / (MU0 / dt + sigma_Hx[j] / 2.0);
        CHX2[j] = 1.0 / (MU0 / dt + sigma_Hx[j] / 2.0) / dy;

        ofs_CHX << j << " " << CHX1[j] << " " << CHX2[j] << " " << sigma_Hx[j] << std::endl;
    }

    for(int i = 0; i < Nx; i++){
        CHY1[i] = (MU0 / dt - sigma_Hy[i] / 2.0) / (MU0 / dt + sigma_Hy[i] / 2.0);
        CHY2[i] = 1.0 / (MU0 / dt + sigma_Hy[i] / 2.0) / dx;              

        ofs_CHY << i << " " << CHY1[i] << " " << CHY2[i] << " " << sigma_Hy[i] << std::endl;
    }

    ofs_CEZX.close();
    ofs_CEZY.close();
    ofs_CHX.close();
    ofs_CHY.close();

    delete [] sigma_Ezx;
    delete [] sigma_Ezy;
    delete [] sigma_Hx;
    delete [] sigma_Hy;
}