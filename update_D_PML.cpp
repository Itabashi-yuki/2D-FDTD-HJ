#include "fdtd2d.h"
#include <iostream>
#include <fstream>

void update_Dx_PML(double ***Dx, double **Hz, double *CDX1, double *CDX2, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 0; i < L; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Dx[NEW][i][j] = CDX1[j] * Dx[OLD][i][j] + CDX2[j] * ( Hz[i][j] - Hz[i][j-1] ); 
        }
    }

    for(int i = Nx - L; i < Nx; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Dx[NEW][i][j] = CDX1[j] * Dx[OLD][i][j] + CDX2[j] * ( Hz[i][j] - Hz[i][j-1] ); 
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = 1; j <= L; j++){
            Dx[NEW][i][j] = CDX1[j] * Dx[OLD][i][j] + CDX2[j] * ( Hz[i][j] - Hz[i][j-1] ); 
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = Ny - L; j <= Ny - 1; j++){
            Dx[NEW][i][j] = CDX1[j] * Dx[OLD][i][j] + CDX2[j] * ( Hz[i][j] - Hz[i][j-1] ); 
        }
    }
}

void update_Dy_PML(double ***Dy, double **Hz, double *CDY1, double *CDY2, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 1; i <= L; i++){
        for(int j = 0; j < Ny; j++){
            Dy[NEW][i][j] = CDY1[i] * Dy[OLD][i][j] - CDY2[i] * ( Hz[i][j] - Hz[i-1][j] );
        }
    }

    for(int i = Nx - L; i <= Nx-1; i++){
        for(int j = 1; j < Ny; j++){
            Dy[NEW][i][j] = CDY1[i] * Dy[OLD][i][j] - CDY2[i] * ( Hz[i][j] - Hz[i-1][j] );
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = 1; j < L; j++){
            Dy[NEW][i][j] = CDY1[i] * Dy[OLD][i][j] - CDY2[i] * ( Hz[i][j] - Hz[i-1][j] );
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = Ny - L; j < Ny; j++){
            Dy[NEW][i][j] = CDY1[i] * Dy[OLD][i][j] - CDY2[i] * ( Hz[i][j] - Hz[i-1][j] );
        }
    }
}

void update_Dz_PML(double ***Dz, double ***Dzx, double ***Dzy, double **Hx, double **Hy, double *CDZX1, double *CDZX2, double *CDZY1, double *CDZY2, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 1; i <= L; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Dzx[NEW][i][j] = CDZX1[i] * Dzx[OLD][i][j] + CDZX2[i] * ( Hy[i][j] - Hy[i-1][j] );
            Dzy[NEW][i][j] = CDZY1[j] * Dzy[OLD][i][j] - CDZY2[j] * ( Hx[i][j] - Hx[i][j-1] );
            Dz[NEW][i][j] = Dzx[NEW][i][j] + Dzy[NEW][i][j];
        }
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Dzx[NEW][i][j] = CDZX1[i] * Dzx[OLD][i][j] + CDZX2[i] * ( Hy[i][j] - Hy[i-1][j] );
            Dzy[NEW][i][j] = CDZY1[j] * Dzy[OLD][i][j] - CDZY2[j] * ( Hx[i][j] - Hx[i][j-1] );
            Dz[NEW][i][j] = Dzx[NEW][i][j] + Dzy[NEW][i][j];
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = 1; j <= L; j++){
            Dzx[NEW][i][j] = CDZX1[i] * Dzx[OLD][i][j] + CDZX2[i] * ( Hy[i][j] - Hy[i-1][j] );
            Dzy[NEW][i][j] = CDZY1[j] * Dzy[OLD][i][j] - CDZY2[j] * ( Hx[i][j] - Hx[i][j-1] );
            Dz[NEW][i][j] = Dzx[NEW][i][j] + Dzy[NEW][i][j];
        } 
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = Ny - L; j <= Ny - 1; j++){
            Dzx[NEW][i][j] = CDZX1[i] * Dzx[OLD][i][j] + CDZX2[i] * ( Hy[i][j] - Hy[i-1][j] );
            Dzy[NEW][i][j] = CDZY1[j] * Dzy[OLD][i][j] - CDZY2[j] * ( Hx[i][j] - Hx[i][j-1] );
            Dz[NEW][i][j] = Dzx[NEW][i][j] + Dzy[NEW][i][j];
        }
    }
}