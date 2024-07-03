#include <iostream>
#include <fstream>
#include "fdtd2d.h"
#include <string>
#include <eigen3/Eigen/Dense>

int main(){

    double **Ex = allocate_2d(Nx, Ny+1, 0.0);
    double **Ey = allocate_2d(Nx+1, Ny, 0.0);
    double **Ez = allocate_2d(Nx+1, Ny+1, 0.0);
    double **Ezx = allocate_2d(Nx+1, Ny+1, 0.0);
    double **Ezy = allocate_2d(Nx+1, Ny+1, 0.0);
    double **Hx = allocate_2d(Nx+1, Ny, 0.0);
    double **Hy = allocate_2d(Nx, Ny+1, 0.0);
    double **Hz = allocate_2d(Nx, Ny, 0.0);
    double **Hzx = allocate_2d(Nx, Ny, 0.0);
    double **Hzy = allocate_2d(Nx, Ny, 0.0);
    double ***Dx = allocate_3d(2, Nx, Ny+1, 0.0);
    double ***Dy = allocate_3d(2, Nx+1, Ny, 0.0);
    double ***Dz = allocate_3d(2, Nx+1, Ny+1, 0.0);
    double ***Dzx = allocate_3d(2, Nx+1, Ny+1, 0.0);
    double ***Dzy = allocate_3d(2, Nx+1, Ny+1, 0.0);
    double ***Jex = allocate_3d(2, Nx, Ny+1, 0.0);
    double ***Jey = allocate_3d(2, Nx+1, Ny, 0.0);
    double ***Jez = allocate_3d(2, Nx+1, Ny+1, 0.0);
    
    double *CEX1 = allocate_1d(Nx, 0.0);
    double *CEX2 = allocate_1d(Nx, 0.0);
    double *CEY1 = allocate_1d(Ny, 0.0);
    double *CEY2 = allocate_1d(Ny, 0.0);
    double *CEZX1 = allocate_1d(Nx, 0.0);
    double *CEZX2 = allocate_1d(Nx, 0.0);
    double *CEZY1 = allocate_1d(Ny, 0.0);
    double *CEZY2 = allocate_1d(Ny, 0.0);
    double *CDX1 = allocate_1d(Nx, 0.0);
    double *CDX2 = allocate_1d(Nx, 0.0);
    double *CDY1 = allocate_1d(Ny, 0.0);
    double *CDY2 = allocate_1d(Ny, 0.0);
    double *CDZX1 = allocate_1d(Nx, 0.0);
    double *CDZX2 = allocate_1d(Nx, 0.0);
    double *CDZY1 = allocate_1d(Ny, 0.0);
    double *CDZY2 = allocate_1d(Ny, 0.0);
    double *CHX1 = allocate_1d(Nx, 0.0);
    double *CHX2 = allocate_1d(Nx, 0.0);
    double *CHY1 = allocate_1d(Ny, 0.0);
    double *CHY2 = allocate_1d(Ny, 0.0);
    double *CHZX1 = allocate_1d(Nx, 0.0);
    double *CHZX2 = allocate_1d(Nx, 0.0);
    double *CHZY1 = allocate_1d(Ny, 0.0);
    double *CHZY2 = allocate_1d(Ny, 0.0);

    Eigen::Matrix3d **S = new Eigen::Matrix3d *[Nx+1];
    Eigen::Matrix3d **B = new Eigen::Matrix3d *[Nx+1];
    Eigen::Matrix3d *S1 = new Eigen::Matrix3d [(Nx+1)*(Ny+1)];
    Eigen::Matrix3d *B1 = new Eigen::Matrix3d [(Nx+1)*(Ny+1)];
    for(int i = 0; i < Nx+1; i++){
        S[i] = S1 + i*(Ny+1);
        B[i] = B1 + i*(Ny+1);
        for(int j = 0; j < Ny+1; j++){
            S[i][j] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
            B[i][j] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        }
    }

    // double Ez_max = -1.0;

    double Ne_exp, nu_exp; /* Neとnuの指数部 10^(Ne_exp) */
    std::cin >> Ne_exp;
    // Ne_exp = 7.0;
    // std::cin >> nu_exp;
    nu_exp = 7.0;

    make_dir();
    output_pal(Ne_exp, nu_exp);
    std::cout << global_dirName << std::endl;

    initialize_Plasma(S, B, Ne_exp, nu_exp);
    initialize_PML(CEX1, CEX2, CEY1, CEY2, CEZX1, CEZX2, CEZY1, CEZY2,
                    CDX1, CDX2, CDY1, CDY2, CDZX1, CDZX2, CDZY1, CDZY2,
                     CHX1, CHX2, CHY1, CHY2, CHZX1, CHZX2, CHZY1, CHZY2);
    // int n0 = cal_obs_n0();
    int n0 = 10;

    // std::ofstream ofs_div_time("./data/" + global_dirName +"/div_time_dt="+ std::to_string(Cdt) + ".dat", std::ios::app);
    std::ofstream ofs_source("./data/" + global_dirName + "/source.dat", std::ios::app);
    for(int n = 1; n < Nt; n++){
        double t = (n - 0.5) * dt;

        if(n % 1000 == 0)
            std::cout << n << " / " << Nt << std::endl;

        update_E(Ex, Ey, Ez, Hx, Hy, Hz, Jex, Jey, Jez, n);
        update_Dx_PML(Dx, Hz, CDX1, CDX2, n);
        update_Dy_PML(Dy, Hz, CDY1, CDY2, n);
        update_Dz_PML(Dz, Dzx, Dzy, Hx, Hy, CDZX1, CDZX2, CDZY1, CDZY2, n);
        // update_Ex_PML(Ex, Hz, CEX1, CEX2);
        // update_Ey_PML(Ey, Hz, CEY1, CEY2);
        // update_Ez_PML(Ez, Ezx, Ezy, Hx, Hy, CEZX1, CEZX2, CEZY1, CEZY2);
        update_Ex_PML(Ex, Dx, Jex, n);
        update_Ey_PML(Ey, Dy, Jey, n);
        update_Ez_PML(Ez, Dz, Jez, n);
        Ez[int(source_x / dx)][int(source_y / dy)] -= dt / EPS0 * Jz(t);

        ofs_source << t << " " << Jz(t) << " " << Ez[int(source_x / dx)][int(source_y / dy)] << std::endl;

        update_H(Hx, Hy, Hz, Ex, Ey, Ez);
        update_Hx_PML(Hx, Ez, CHX1, CHX2);
        update_Hy_PML(Hy, Ez, CHY1, CHY2);
        update_Hz_PML(Hz, Hzx, Hzy, Ex, Ey, CHZX1, CHZX2, CHZY1, CHZY2);

        update_J(Jex, Jey, Jez, Ex, Ey, Ez, S, B, n);
        
        // output_Ez(Ex, Ey, Ez, n, n0);
        output_obs(Ex, Ey, Ez, Hx, Hy, Hz, Jex, Jey, Jez, n);

        // if(std::abs(Ez[int(obs_x / dx)][int(obs_y / dy)]) > 10000){
            // ofs_div_time << return_Ne(Ne_alt) << " " << cal_nu(nu_alt) << " " << n * dt << std::endl;
            // exit(0);
        // }
        // judge_Ez_div(Ez, n, n0, &Ez_max, Ne_alt, nu_alt);
    }

    // ofs_div_time.close();


    free_memory3d(Jex, 2, Nx);
    free_memory3d(Jey, 2, Nx+1);
    free_memory3d(Jez, 2, Nx+1);
    free_memory2d(Ex, Nx);
    free_memory2d(Ey, Nx+1);
    free_memory2d(Ez, Nx+1);
    free_memory2d(Ezx, Nx+1);
    free_memory2d(Ezy, Nx+1);
    free_memory3d(Dx, 2, Nx);
    free_memory3d(Dy, 2, Ny+1);
    free_memory3d(Dz, 2, Nx+1);
    free_memory3d(Dzx, 2, Nx+1);
    free_memory3d(Dzy, 2, Nx+1);
    free_memory2d(Hx, Nx+1);
    free_memory2d(Hy, Nx);
    free_memory2d(Hz, Nx);
    free_memory2d(Hzx, Nx);
    free_memory2d(Hzy, Nx);

    delete [] CEX1;
    delete [] CEX2;
    delete [] CEY1;
    delete [] CEY2;
    delete [] CEZX1;
    delete [] CEZX2;
    delete [] CEZY1;
    delete [] CEZY2;
    delete [] CDX1;
    delete [] CDX2;
    delete [] CDY1;
    delete [] CDY2;
    delete [] CDZX1;
    delete [] CDZX2;
    delete [] CDZY1;
    delete [] CDZY2;
    delete [] CHX1;
    delete [] CHX2;
    delete [] CHY1;
    delete [] CHY2;
    delete [] CHZX1;
    delete [] CHZX2;
    delete [] CHZY1;
    delete [] CHZY2;
}