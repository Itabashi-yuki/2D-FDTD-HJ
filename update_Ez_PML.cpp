#include "fdtd2d.h"
#include <math.h>
#include <iostream>
#include <fstream>

void update_E_PML(double **Ez, double **Ezx, double **Ezy, double **Hx, double **Hy, double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2){
    // PML_domain_E domain;

    // for(int pml_idx = 0; pml_idx < 4; pml_idx++){
    //     for(int i = domain.idx[pml_idx].i1; i <= domain.idx[pml_idx].i2; i++){
    //         for(int j = domain.idx[pml_idx].j1; j <= domain.idx[pml_idx].j2; j++){
    //             Ezx[i][j] = CEZX1[i] * Ezx[i][j] + CEZX2[i] * (Hy[i][j] - Hy[i-1][j]);
    //             Ezy[i][j] = CEZY1[i] * Ezy[i][j] - CEZY2[i] * (Hx[i][j] - Hx[i][j-1]);
    //             Ez[i][j] = Ezx[i][j] + Ezy[i][j];
    //         }
    //     }
    // }

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