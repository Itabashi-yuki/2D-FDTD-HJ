#include "fdtd2d.h"

void update_E(double **Ex, double **Ey, double **Ez, double **Hx, double **Hy, double **Hz,
                 double ***Jex, double ***Jey, double ***Jez, int n){
    int NEW = n % 2;
    int OLD = (n+1) % 2;

    for(int i = L; i < Nx - L; i++){
        for(int j = L + 1; j <= Ny - L - 1; j++){
            // Ex[i][j] = Ex[i][j] + dt / EPS0 / dy * ( Hz[i][j] - Hz[i][j-1] ) - dt / EPS0 * Jex[NEW][i][j];
            Ex[i][j] = Ex[i][j] + dt / EPS0 / dy * ( Hz[i][j] - Hz[i][j-1] ) - dt / EPS0 * Jex[OLD][i][j];
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = L; j < Ny - L; j++){
            // Ey[i][j] = Ey[i][j] - dt / EPS0 / dx * ( Hz[i][j] - Hz[i-1][j] ) - dt / EPS0 * Jey[NEW][i][j];
            Ey[i][j] = Ey[i][j] - dt / EPS0 / dx * ( Hz[i][j] - Hz[i-1][j] ) - dt / EPS0 * Jey[OLD][i][j];
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = L + 1; j <= Ny - L - 1; j++){
            // Ez[i][j] = Ez[i][j] + dt / EPS0 / dx * ( Hy[i][j] - Hy[i-1][j] ) - dt / EPS0 / dy * ( Hx[i][j] - Hx[i][j-1] )
            //  - dt / EPS0 * Jez[NEW][i][j];
            Ez[i][j] = Ez[i][j] + dt / EPS0 / dx * ( Hy[i][j] - Hy[i-1][j] ) - dt / EPS0 / dy * ( Hx[i][j] - Hx[i][j-1] )
             - dt / EPS0 * Jez[OLD][i][j];
        }
    }


    // Ez[int(source_x / dx)][int(source_y / dy)] = Ez[int(source_x / dx)][int(source_y / dy)] - dt / EPS0 * Jz[n-1];
}