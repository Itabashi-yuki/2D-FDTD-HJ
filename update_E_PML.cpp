#include "fdtd2d.h"
#include <math.h>
#include <iostream>
#include <fstream>

/* Ex の 更新式 */
void update_Ex_PML(double **Ex, double **Hz, double *CEX1, double *CEX2){
    for(int i = 0; i < L; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Ex[i][j] = CEX1[j] * Ex[i][j] + CEX2[j] * (Hz[i][j] - Hz[i][j-1]);
        }
    }

    for(int i = Nx - L; i < Nx; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Ex[i][j] = CEX1[j] * Ex[i][j] + CEX2[j] * (Hz[i][j] - Hz[i][j-1]);
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = 1; j <= L; j++){
            Ex[i][j] = CEX1[j] * Ex[i][j] + CEX2[j] * (Hz[i][j] - Hz[i][j-1]);
        }
    }

    for(int i = L; i < Nx - L; i++){
        for(int j = Ny - L; j <= Ny - 1; j++){
            Ex[i][j] = CEX1[j] * Ex[i][j] + CEX2[j] * (Hz[i][j] - Hz[i][j-1]);
        }
    }
}

/* Ey の 更新式 */
void update_Ey_PML(double **Ey, double **Hz, double *CEY1, double *CEY2){
    for(int i = 1; i <= L; i++){
        for(int j = 0; j < Ny; j++){
            Ey[i][j] = CEY1[i] * Ey[i][j] - CEY2[i] * (Hz[i][j] - Hz[i-1][j]);
        }
    }

    for(int i = Nx - L; i <= Nx-1; i++){
        for(int j = 1; j < Ny; j++){
            Ey[i][j] = CEY1[i] * Ey[i][j] - CEY2[i] * (Hz[i][j] - Hz[i-1][j]);
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = 1; j < L; j++){
            Ey[i][j] = CEY1[i] * Ey[i][j] - CEY2[i] * (Hz[i][j] - Hz[i-1][j]);
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = Ny - L; j < Ny; j++){
            Ey[i][j] = CEY1[i] * Ey[i][j] - CEY2[i] * (Hz[i][j] - Hz[i-1][j]);
        }
    }
}

/* Ez の 更新式 */
void update_Ez_PML(double **Ez, double **Ezx, double **Ezy, double **Hx, double **Hy,
                    double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2){

    /* Ez の 更新式 */
    for(int i = 1; i <= L; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Ezx[i][j] = CEZX1[i] * Ezx[i][j] + CEZX2[i] * (Hy[i][j] - Hy[i-1][j]);
            Ezy[i][j] = CEZY1[j] * Ezy[i][j] - CEZY2[j] * (Hx[i][j] - Hx[i][j-1]);
            Ez[i][j] = Ezx[i][j] + Ezy[i][j];
        }
    }

    for(int i = Nx - L; i <= Nx - 1; i++){
        for(int j = 1; j <= Ny - 1; j++){
            Ezx[i][j] = CEZX1[i] * Ezx[i][j] + CEZX2[i] * (Hy[i][j] - Hy[i-1][j]);
            Ezy[i][j] = CEZY1[j] * Ezy[i][j] - CEZY2[j] * (Hx[i][j] - Hx[i][j-1]);
            Ez[i][j] = Ezx[i][j] + Ezy[i][j];
        }
    }

    for(int i = L + 1; i <= Nx - L - 1; i++){
        for(int j = 1; j <= L; j++){
            Ezx[i][j] = CEZX1[i] * Ezx[i][j] + CEZX2[i] * (Hy[i][j] - Hy[i-1][j]);
            Ezy[i][j] = CEZY1[j] * Ezy[i][j] - CEZY2[j] * (Hx[i][j] - Hx[i][j-1]);
            Ez[i][j] = Ezx[i][j] + Ezy[i][j];
        }
        for(int j = Ny - L; j <= Ny - 1; j++){
            Ezx[i][j] = CEZX1[i] * Ezx[i][j] + CEZX2[i] * (Hy[i][j] - Hy[i-1][j]);
            Ezy[i][j] = CEZY1[j] * Ezy[i][j] - CEZY2[j] * (Hx[i][j] - Hx[i][j-1]);
            Ez[i][j] = Ezx[i][j] + Ezy[i][j];
        }
        
    }
}