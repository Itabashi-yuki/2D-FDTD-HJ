#include "fdtd2d.h"
#include <math.h>
#include <iostream>
#include <fstream>

void update_H_PML(double **Hx, double **Hy, double **Ez, double *CHX1, double *CHX2, double *CHY1, double *CHY2){
    // PML_domain_Hx domain_Hx;
    // PML_domain_Hy domain_Hy;
    // for(int pml_idx = 0; pml_idx < 4; pml_idx++){
    //     for(int i = domain_Hx.idx[pml_idx].i1; i <= domain_Hx.idx[pml_idx].i2; i++){
    //         for(int j = domain_Hx.idx[pml_idx].j1; j <= domain_Hx.idx[pml_idx].j2; j++){
    //             Hx[i][j] = CHX1[i] * Hx[i][j] - CHX2[i] * (Ez[i][j+1] - Ez[i][j]);
    //         }
    //     }
    // }

    // for(int pml_idx = 0; pml_idx < 4; pml_idx++){
    //     for(int i = domain_Hy.idx[pml_idx].i1; i <= domain_Hy.idx[pml_idx].i2; i++){
    //         for(int j = domain_Hy.idx[pml_idx].j1; j <= domain_Hy.idx[pml_idx].j2; j++){
    //             Hy[i][j] = CHY1[i] * Hy[i][j] + CHY2[i] * (Ez[i+1][j] - Ez[i][j]);
    //         }
    //     }
    // }

    for(int i = 1; i <= L; i++){
        for(int j = 0; j < Ny; j++){
            Hx[i][j] = CHX1[j] * Hx[i][j] - CHX2[j] * (Ez[i][j+1] - Ez[i][j]);
        }
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        for(int j = 0; j < Ny; j++){
            Hx[i][j] = CHX1[j] * Hx[i][j] - CHX2[j] * (Ez[i][j+1] - Ez[i][j]);
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = 0; j < L; j++){
            Hx[i][j] = CHX1[j] * Hx[i][j] - CHX2[j] * (Ez[i][j+1] - Ez[i][j]);
        } 
        for(int j = Ny - L; j <= Ny - 1; j++){
            Hx[i][j] = CHX1[j] * Hx[i][j] - CHX2[j] * (Ez[i][j+1] - Ez[i][j]);
        }
    }



    for(int i = 0; i < L; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Hy[i][j] = CHY1[i] * Hy[i][j] + CHY2[i] * (Ez[i+1][j] - Ez[i][j]);
        }
    }

    for(int i = Nx - L; i < Nx; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Hy[i][j] = CHY1[i] * Hy[i][j] + CHY2[i] * (Ez[i+1][j] - Ez[i][j]);
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = 1; j <= L; j++){
            Hy[i][j] = CHY1[i] * Hy[i][j] + CHY2[i] * (Ez[i+1][j] - Ez[i][j]);
        }
        for(int j = Ny - L; j <= Ny - 1; j++){
            Hy[i][j] = CHY1[i] * Hy[i][j] + CHY2[i] * (Ez[i+1][j] - Ez[i][j]);
        }
    }
}