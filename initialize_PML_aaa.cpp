#include "fdtd2d.h"
#include <math.h>
#include <iostream>
#include <fstream>

void initialize_PML(double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2, double *CHX1, double *CHX2, double *CHY1, double *CHY2){
    double sigma_max = - (M + 1.0) * Y / 2.0 / L / dx * std::log(R);
    double *sigma_Ex = new double[Nx];
    double *sigma_Ey = new double[Ny];
    double *sigma_Hx = new double[Nx];
    double *sigma_Hy = new double[Nx];

    PML_domain_E domain;    
    PML_domain_Hx domain_Hx;
    PML_domain_Hy domain_Hy;

    // for(int pml_idx = 0; pml_idx < 4; pml_idx++){
    //     std::cout << "i: " <<domain_Hx.idx[pml_idx].i1 << " ~ " << domain_Hx.idx[pml_idx].i2 << std::endl;
    //     std::cout << "j: " <<domain_Hx.idx[pml_idx].j1 << " ~ " << domain_Hx.idx[pml_idx].j2 << std::endl;
    //     std::cout << std::endl;
    // }
    

    // exit(0);

    

    for(int pml_idx = 0; pml_idx < 3; pml_idx++){
        for(int i = domain.idx[pml_idx].i1; i <= domain.idx[pml_idx].i2; i++){
            if(pml_idx == 0){
                sigma_Ex[i] = sigma_max * std::pow(((L*dx - i*dx) / L / dx), M);
            } else if(pml_idx == 1) {
                sigma_Ex[i] = sigma_max * std::pow(((i*dx - (Nx - L) * dx) / L / dx), M);
            } else if(pml_idx == 2){
                sigma_Ex[i] = 0.0;
            }
            CEZX1[i] = (EPS0 / dt - sigma_Ex[i] / 2.0) / (EPS0 / dt + sigma_Ex[i] / 2.0);
            CEZX2[i] = 1.0 / (EPS0 / dt + sigma_Ex[i] / 2.0) / dx;
        }
    }

    for(int pml_idx = 1; pml_idx < 4; pml_idx++){
        for(int i = domain.idx[pml_idx].j1; i <= domain.idx[pml_idx].j2; i++){
            if(pml_idx == 1){
                sigma_Ey[i] = 0.0;
            } else if(pml_idx == 2){
                sigma_Ey[i] =sigma_max * std::pow(((L*dy - i*dy) / L / dy), M);
            } else if(pml_idx == 3){
                sigma_Ey[i] = sigma_max * std::pow(((i*dy - (Ny - L) * dy) / L / dy), M);
            }
            CEZY1[i] = (EPS0 / dt - sigma_Ey[i] / 2.0) / (EPS0 / dt + sigma_Ey[i] / 2.0);
            CEZY2[i] = 1.0 / (EPS0 / dt + sigma_Ey[i] / 2.0) / dy;
        }
    }

    std::ofstream ofs_CEZX("CEZX.dat");
    std::ofstream ofs_CEZY("CEZY.dat");

    for(int i = 0; i < Nx; i++){
        ofs_CEZX << i << " " << CEZX1[i] << " " << CEZX2[i] << " " << sigma_Ex[i] << std::endl;
        ofs_CEZY << i << " " << CEZY1[i] << " " << CEZY2[i] << " " << sigma_Ey[i] << std::endl;
    }

    for(int pml_idx = 1; pml_idx < 4; pml_idx++){
        for(int i = domain_Hx.idx[pml_idx].j1; i <= domain_Hx.idx[pml_idx].j2; i++){
            if(pml_idx == 1){
                sigma_Hy[i] = 0.0;
            } else if(pml_idx == 2) {
                sigma_Hy[i] = Z * Z * sigma_max * std::pow(((L*dy - (i + 0.5)*dy) / L / dy), M);
            } else if(pml_idx == 3){
                sigma_Hy[i] = Z * Z * sigma_max * std::pow((((i + 0.5)*dy - (Nx - L) * dy) / L / dy), M);
            }
            CHX1[i] = (MU0 / dt - sigma_Hy[i] / 2.0) / (MU0 / dt + sigma_Hy[i] / 2.0);
            CHX2[i] = 1.0 / (MU0 / dt + sigma_Hy[i] / 2.0) / dy;
        }
    }

    for(int pml_idx = 0; pml_idx < 3; pml_idx++){
        for(int i = domain_Hy.idx[pml_idx].i1; i <= domain_Hy.idx[pml_idx].i2; i++){
            if(pml_idx == 0){
                sigma_Hx[i] = Z * Z * sigma_max * std::pow(((L*dx - (i + 0.5)*dx) / L / dx), M);
            } else if(pml_idx == 1){
                sigma_Hx[i] = Z * Z * sigma_max * std::pow((((i + 0.5)*dx - (Nx - L) * dx) / L / dx), M);
            } else if(pml_idx == 2){
                sigma_Hx[i] = 0.0;
            }
            CHY1[i] = (MU0 / dt - sigma_Hx[i] / 2.0) / (MU0 / dt + sigma_Hx[i] / 2.0);
            CHY2[i] = 1.0 / (MU0 / dt + sigma_Hx[i] / 2.0) / dx;
        }
    }

    std::ofstream ofs_CHX("CHX.dat");
    std::ofstream ofs_CHY("CHY.dat");

    for(int i = 0; i < Nx; i++){
        ofs_CHX << i << " " << CHX1[i] << " " << CHX2[i] << " " << sigma_Hy[i] << std::endl;
        ofs_CHY << i << " " << CHY1[i] << " " << CHY2[i] << " " << sigma_Hx[i] << std::endl;
    }

    delete [] sigma_Ex;
    delete [] sigma_Ey;
    delete [] sigma_Hx;
    delete [] sigma_Hy;
}