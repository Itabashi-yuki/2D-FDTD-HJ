#include "fdtd2d.h"
#include <fstream>
#include <iostream>

void initialize_PML(double *CEX1, double *CEX2, double *CEY1, double *CEY2,
                    double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2,
                    double *CDX1, double *CDX2, double *CDY1, double *CDY2,
                    double *CDZX1, double *CDZX2, double *CDZY1, double *CDZY2,  
                    double *CHX1, double *CHX2, double *CHY1, double *CHY2,
                    double *CHZX1, double *CHZX2, double *CHZY1, double *CHZY2
                    ){
    // double sigma_max_x = - (M + 1.0) * Y / 2.0 / L / dx * std::log(R);
    // double sigma_max_y = - (M + 1.0) * Y / 2.0 / L / dy * std::log(R);
    double sigma_max_x = - (M + 1.0) * C0 / 2.0 / L / dx * std::log(R); /* Complex coordinate stretching によるsigma_max_x */
    double sigma_max_y = - (M + 1.0) * C0 / 2.0 / L / dy * std::log(R); /* Complex coordinate stretching によるsigma_max_y */
    double *sigma_Ex = allocate_1d(Nx, 0.0);
    double *sigma_Ey = allocate_1d(Ny, 0.0);
    double *sigma_Ezx = allocate_1d(Ny, 0.0);
    double *sigma_Ezy = allocate_1d(Nx, 0.0);
    double *sigma_Dx = allocate_1d(Nx, 0.0);
    double *sigma_Dy = allocate_1d(Ny, 0.0);
    double *sigma_Dzx = allocate_1d(Nx, 0.0);
    double *sigma_Dzy = allocate_1d(Ny, 0.0);
    double *sigma_Hx = allocate_1d(Ny, 0.0);
    double *sigma_Hy = allocate_1d(Nx, 0.0);
    double *sigma_Hzx = allocate_1d(Nx, 0.0);
    double *sigma_Hzy = allocate_1d(Ny, 0.0);

    /* σ の設定 */

    for(int j = 1; j <= L; j++){
        sigma_Ex[j] = sigma_max_y * std::pow(((L*dy - j*dy) / L / dy), M);
    }

    for(int i = 1; i <= L; i++){
        sigma_Ey[i] = sigma_max_x * std::pow(((L*dx - i*dx) / L / dx), M);
    }

    for(int i = 1; i <= L; i++){
        sigma_Ezx[i] = sigma_max_x * std::pow(((L*dx - i*dx) / L / dx), M);
    }

    for(int j = 1; j <= L; j++){
        sigma_Ezy[j] = sigma_max_y * std::pow(((L*dy - j*dy) / L / dy), M);
    }

    for(int j = 1; j <= L; j++){
        sigma_Dx[j] = sigma_max_y * std::pow(((L*dy - j*dy) / L / dy), M);
    }

    for(int i = 1; i <= L; i++){
        sigma_Dy[i] = sigma_max_x * std::pow(((L*dx - i*dx) / L / dx), M);
    }

    for(int i = 1; i <= L; i++){
        sigma_Dzx[i] = sigma_max_x * std::pow(((L*dx - i*dx) / L / dx), M);
    }

    for(int j = 1; j <= L; j++){
        sigma_Dzy[j] = sigma_max_y * std::pow(((L*dy - j*dy) / L / dy), M);
    }

    for(int j = 0; j < L; j++){
        sigma_Hx[j] = sigma_max_y * std::pow(((L*dy - (j + 0.5) * dy) / L / dy), M);
    }

    for(int i = 0; i < L; i++){
        sigma_Hy[i] = sigma_max_x * std::pow(((L*dx - (i + 0.5)*dx) / L / dx), M);
    }

    for(int i = 0; i < L; i++){
        sigma_Hzx[i] = sigma_max_x * std::pow(((L*dx - (i + 0.5)*dx) / L / dx), M);
    }

    for(int j = 0; j < L; j++){
        sigma_Hzy[j] = sigma_max_y * std::pow(((L*dy - (j + 0.5) * dy) / L / dy), M);
    }

    for(int j = Ny - L; j <= Ny - 1; j++){
        sigma_Ex[j] = sigma_max_y * std::pow(((j*dy -(Ny - L) * dy) / L / dy), M);
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        sigma_Ey[i] = sigma_max_x * std::pow(((i*dx -(Nx - L) * dx) / L / dx), M);
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        sigma_Ezx[i] = sigma_max_x * std::pow(((i*dx - (Nx - L) * dx) / L / dx), M);
    }

    for(int j = Ny - L; j <= Ny - 1; j++){
        sigma_Ezy[j] = sigma_max_y * std::pow(((j*dy - (Ny - L) * dy) / L / dy), M);
    }

    for(int j = Ny - L; j <= Ny - 1; j++){
        sigma_Dx[j] = sigma_max_y * std::pow(((j*dy -(Ny - L) * dy) / L / dy), M);
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        sigma_Dy[i] = sigma_max_x * std::pow(((i*dx -(Nx - L) * dx) / L / dx), M);
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        sigma_Dzx[i] = sigma_max_x * std::pow(((i*dx - (Nx - L) * dx) / L / dx), M);
    }

    for(int j = Ny - L; j <= Ny - 1; j++){
        sigma_Dzy[j] = sigma_max_y * std::pow(((j*dy - (Ny - L) * dy) / L / dy), M);
    }

    for(int j = Ny - L; j < Ny; j++){
        sigma_Hx[j] = sigma_max_y * std::pow((((j + 0.5) * dy - (Ny - L) * dy) / L / dy), M);
    }

    for(int i = Nx - L; i < Nx; i++){
        sigma_Hy[i] = sigma_max_x * std::pow((((i + 0.5) * dx - (Nx - L) * dx) / L / dx), M);
    }

    for(int i = Nx - L; i < Nx; i++){
        sigma_Hzx[i] = sigma_max_x * std::pow((((i + 0.5) * dx - (Nx - L) * dx) / L / dx), M);
    }

    for(int j = Ny - L; j < Ny; j++){
        sigma_Hzy[j] = sigma_max_y * std::pow((((j + 0.5) * dy - (Ny - L) * dy) / L / dy), M);
    }

    std::ofstream ofs_sigma("./data/" + global_dirName + "/sigma.dat");
    for(int i = 0; i < Nx; i++){
        ofs_sigma << sigma_Dzx[i] << " " << sigma_Hzx[i] << std::endl;
    }
    // exit(0);

    for(int j = 1; j <= Ny - 1; j++){
        CEX1[j] = (1.0 / dt - sigma_Ex[j] / 2.0) / (1.0 / dt + sigma_Ex[j] / 2.0);
        CEX2[j] = 1.0 / (1.0 / dt + sigma_Ex[j] / 2.0) / EPS0 / dy;
    }

    for(int i = 1; i <= Nx - 1; i++){
        CEY1[i] = (1.0 / dt - sigma_Ey[i] / 2.0) / (1.0 / dt + sigma_Ey[i] / 2.0);
        CEY2[i] = 1.0 / (1.0 / dt + sigma_Ey[i] / 2.0) / EPS0 / dx;
    }

    for(int i = 1; i <= Nx - 1; i++){
        CEZX1[i] = (1.0 / dt - sigma_Ezx[i] / 2.0) / (1.0 / dt + sigma_Ezx[i] / 2.0);
        CEZX2[i] = 1.0 / (1.0 / dt + sigma_Ezx[i] / 2.0) / EPS0 / dx;
    }

    for(int j = 1; j <= Ny - 1; j++){
        CEZY1[j] = (1.0 / dt - sigma_Ezy[j] / 2.0) / (1.0 / dt + sigma_Ezy[j] / 2.0);
        CEZY2[j] = 1.0 / (1.0 / dt + sigma_Ezy[j] / 2.0) / EPS0 / dy;
    }

    for(int j = 1; j <= Ny - 1; j++){
        CDX1[j] = (1.0 / dt - sigma_Dx[j] / 2.0) / (1.0 / dt + sigma_Dx[j] / 2.0);
        CDX2[j] = 1.0 / (1.0 / dt + sigma_Dx[j] / 2.0) / dy;
    }

    for(int i = 1; i <= Nx - 1; i++){
        CDY1[i] = (1.0 / dt - sigma_Dy[i] / 2.0) / (1.0 / dt + sigma_Dy[i] / 2.0);
        CDY2[i] = 1.0 / (1.0 / dt + sigma_Dy[i] / 2.0) / dx;
    }

    for(int i = 1; i <= Nx - 1; i++){
        CDZX1[i] = (1.0 / dt - sigma_Dzx[i] / 2.0) / (1.0 / dt + sigma_Dzx[i] / 2.0);
        CDZX2[i] = 1.0 / (1.0 / dt + sigma_Dzx[i] / 2.0) / dx;
    }

    for(int j = 1; j <= Ny - 1; j++){
        CDZY1[j] = (1.0 / dt - sigma_Dzy[j] / 2.0) / (1.0 / dt + sigma_Dzy[j] / 2.0);
        CDZY2[j] = 1.0 / (1.0 / dt + sigma_Dzy[j] / 2.0) / dy;
    }

    for(int j = 0; j < Ny; j++){
        CHX1[j] = (1.0 / dt - sigma_Hx[j] / 2.0) / (1.0 / dt + sigma_Hx[j] / 2.0);
        CHX2[j] = 1.0 / (1.0 / dt + sigma_Hx[j] / 2.0) / MU0 / dy;
    }

    for(int i = 0; i < Nx; i++){
        CHY1[i] = (1.0 / dt - sigma_Hy[i] / 2.0) / (1.0 / dt + sigma_Hy[i] / 2.0);
        CHY2[i] = 1.0 / (1.0 / dt + sigma_Hy[i] / 2.0) / MU0 / dx;   
    }

    for(int i = 0; i < Nx; i++){
        CHZX1[i] = (1.0 / dt - sigma_Hzx[i] / 2.0) / (1.0 / dt + sigma_Hzx[i] / 2.0);
        CHZX2[i] = 1.0 / (1.0 / dt + sigma_Hzx[i] / 2.0) / MU0 / dx;
    }

    for(int j = 0; j < Ny; j++){
        CHZY1[j] = (1.0 / dt - sigma_Hzy[j] / 2.0) / (1.0 / dt + sigma_Hzy[j] / 2.0);
        CHZY2[j] = 1.0 / (1.0 / dt + sigma_Hzy[j] / 2.0) / MU0 / dy;
    }

    std::ofstream ofs_C("./data/" + global_dirName + "/C.dat");
    for(int i = 0; i < Nx; i++){
        ofs_C << i << " " << CEX1[i] << " " << CEX2[i] << " " << CEY1[i] << " " << CEY2[i] << " "
                << CEZX1[i] << " " << CEZX2[i] << " " << CEZY1[i] << " " << CEZY2[i] << " " 
                << CDX1[i] << " " << CDX2[i] << " " << CDY1[i] << " " << CDY2[i] << " " 
                << CDZX1[i] << " " << CDZX2[i] << " " << CDZY1[i] << " " << CDZY2[i] << " "
                << CHX1[i] << " " << CHX2[i] <<" " << CHY1[i] << " " << CHY2[i] << " " 
                << CHZX1[i] << " " << CHZX2[i] << " " << CHZY1[i] << " " << CHZY2[i] << std::endl;
    }

    // exit(0);

    delete [] sigma_Ex;
    delete [] sigma_Ey;
    delete [] sigma_Ezx;
    delete [] sigma_Ezy;
    delete [] sigma_Dx;
    delete [] sigma_Dy;
    delete [] sigma_Dzx;
    delete [] sigma_Dzy;
    delete [] sigma_Hx;
    delete [] sigma_Hy;
    delete [] sigma_Hzx;
    delete [] sigma_Hzy;
}