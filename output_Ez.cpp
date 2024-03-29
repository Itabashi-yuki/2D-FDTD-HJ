#include <iostream>
#include <fstream>
#include <string>
#include "fdtd2d.h"


void output_Ez(double **Ex, double **Ey, double **Ez, int n, int n0){
    if(n == 1 || n % n0 == 0 || n == Nt - 1){
        std::ofstream ofs_Ez("./data/" + global_dirName + "/E_" + std::to_string(n) + ".dat");

        ofs_Ez << "#" << n0 * dt << "[s]ごと出力, 以下は" << n*dt << "[s]における結果" <<  std::endl;

        for(int i = 0; i < Nx; i+=2){
            for(int j = 0; j < Ny; j+=2){
                // ofs_Ez << i * dx * 1e-3 << " " << j * dy * 1e-3 << " " << Ez[i][j] << std::endl;
                ofs_Ez << i * dx * 1e-3 << " " << j * dy * 1e-3 << " " 
                << Ex[i][j] << " " << Ey[i][j] << " " << Ez[i][j] << " " 
                << std::sqrt(Ex[i][j]*Ex[i][j] + Ey[i][j]*Ey[i][j] + Ez[i][j]*Ez[i][j]) << std::endl;
            }
            ofs_Ez << std::endl;
        }

    ofs_Ez.close();

    }
}