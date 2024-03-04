#include "fdtd2d.h"

void update_E(double **Ez, double **Hx, double **Hy, double *Jz, int n){
    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = L + 1; j <= Ny - L - 1; j++){
            Ez[i][j] = Ez[i][j] + dt / EPS0 / dx * ( Hy[i][j] - Hy[i-1][j] ) - dt / EPS0 / dy * ( Hx[i][j] - Hx[i][j-1] );
        }
    }
    Ez[Nx/2][Ny/2] = Ez[Nx/2][Ny/2] - dt / EPS0 * Jz[n-1];
}