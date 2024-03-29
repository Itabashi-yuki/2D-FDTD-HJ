#include "fdtd2d.h"
#include <math.h>
#include <iostream>
#include <fstream>

void update_Hx_PML(double **Hx, double **Ez, double *CHX1, double *CHX2){
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
        for(int j = Ny - L; j < Ny; j++){
            Hx[i][j] = CHX1[j] * Hx[i][j] - CHX2[j] * (Ez[i][j+1] - Ez[i][j]);
        }
    }
}

void update_Hy_PML(double **Hy, double **Ez, double *CHY1, double *CHY2){
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

void update_Hz_PML(double **Hz, double **Hzx, double **Hzy, 
                 double **Ex, double **Ey,double *CHZX1, double *CHZX2,
                 double *CHZY1, double *CHZY2){
    for(int i = 0; i < L; i++){
        for(int j = 0; j < Ny; j++){
            Hzx[i][j] = CHZX1[i] * Hzx[i][j] - CHZX2[i] * (Ey[i+1][j] - Ey[i][j]);
            Hzy[i][j] = CHZY1[j] * Hzy[i][j] + CHZY2[j] * (Ex[i][j+1] - Ex[i][j]);
            Hz[i][j] = Hzx[i][j] + Hzy[i][j];
        }
    }

    for(int i = Nx - L; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            Hzx[i][j] = CHZX1[i] * Hzx[i][j] - CHZX2[i] * (Ey[i+1][j] - Ey[i][j]);
            Hzy[i][j] = CHZY1[j] * Hzy[i][j] + CHZY2[j] * (Ex[i][j+1] - Ex[i][j]);
            Hz[i][j] = Hzx[i][j] + Hzy[i][j];
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = 0; j < L; j++){
            Hzx[i][j] = CHZX1[i] * Hzx[i][j] - CHZX2[i] * (Ey[i+1][j] - Ey[i][j]);
            Hzy[i][j] = CHZY1[j] * Hzy[i][j] + CHZY2[j] * (Ex[i][j+1] - Ex[i][j]);
            Hz[i][j] = Hzx[i][j] + Hzy[i][j];
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = Ny - L; j < Ny; j++){
            Hzx[i][j] = CHZX1[i] * Hzx[i][j] - CHZX2[i] * (Ey[i+1][j] - Ey[i][j]);
            Hzy[i][j] = CHZY1[j] * Hzy[i][j] + CHZY2[j] * (Ex[i][j+1] - Ex[i][j]);
            Hz[i][j] = Hzx[i][j] + Hzy[i][j];
        }
    }
}