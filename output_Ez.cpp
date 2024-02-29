#include <iostream>
#include <fstream>
#include <string>
#include "fdtd2d.h"

void output_Ez(double **Ez, int n){
    std::ofstream ofs_Ez("test/Ez_" + std::to_string(n) + ".dat");

    for(int i = 0; i < Nx + 1; i++){
        for(int j = 0; j < Ny + 1; j++){
            ofs_Ez << i * dx * 1e-3 << " " << j * dy * 1e-3 << " " << Ez[i][j] << std::endl;
        }
        ofs_Ez << std::endl;
    }
    ofs_Ez.close();
}