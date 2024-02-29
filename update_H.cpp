#include "fdtd2d.h"

void update_H(double **Hx, double **Hy, double **Ez){
    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = L; j < Ny - L; j++){
            Hx[i][j] = Hx[i][j] - dt / MU0 / dy * (Ez[i][j+1] - Ez[i][j]);
       }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = L + 1; j <= Ny - L - 1; j++){
            Hy[i][j] = Hy[i][j] + dt / MU0 / dx * (Ez[i+1][j] - Ez[i][j]);
        }
    }
}