#include <iostream>
#include <fstream>
#include "fdtd2d.h"

int main(){
    double **Ez = new double *[Nx+1];
    double **Ezx = new double *[Nx+1];
    double **Ezy = new double *[Nx+1];
    double **Hx = new double *[Nx+1];
    double **Hy = new double *[Nx+1];
    double *CEZX1 = new double[Nx];
    double *CEZX2 = new double[Nx];
    double *CEZY1 = new double[Ny];
    double *CEZY2 = new double[Ny];
    double *CHX1 = new double[Nx];
    double *CHX2 = new double[Nx];
    double *CHY1 = new double[Ny];
    double *CHY2 = new double[Ny];
    for(int i = 0; i < Nx+1; i++){
        Ez[i] = new double [Ny+1];
        Ezx[i] = new double [Ny+1];
        Ezy[i] = new double [Ny+1];
        Hx[i] = new double [Ny+1];
        Hy[i] = new double [Ny+1];
    }
    double *Jz = new double [Nt];
    
    initialize_PML(CEZX1, CEZX2, CEZY1, CEZY2, CHX1, CHX2, CHY1, CHY2);

    std::ofstream ofs("obs2.dat");
    for(int n = 1; n < Nt; n++){
        if(n % 100 == 0)
            std::cout << n << " / " << Nt << std::endl;
        current_source(Jz, n);
        update_E(Ez, Hx, Hy, Jz, n);
        update_E_PML(Ez, Ezx, Ezy, Hx, Hy, CEZX1, CEZX2, CEZY1, CEZY2);
        update_H(Hx, Hy, Ez);
        update_H_PML(Hx, Hy, Ez, CHX1, CHX2, CHY1, CHY2);
        output_Ez(Ez, n);
        ofs << n * dt << " " << Ez[75][50] << std::endl;
    }

    ofs.close();

    for(int i = 0; i < Nx+1; i++){
        delete [] Ez[i];
        delete [] Ezx[i];
        delete [] Ezy[i];
        delete [] Hx[i];
        delete [] Hy[i];
    }
    delete [] Jz;
    delete [] CEZX1;
    delete [] CEZX2;
    delete [] CEZY1;
    delete [] CEZY2;
    delete [] CHX1;
    delete [] CHX2;
    delete [] CHY1;
    delete [] CHY2;
}