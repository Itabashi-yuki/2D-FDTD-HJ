#include "fdtd2d.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>

void update_J(double ***Jex, double ***Jey, double ***Jez, double **Ex, double **Ey, double **Ez, Eigen::Matrix3d **S, Eigen::Matrix3d **B, int n){
    int NEW = n % 2;
    int OLD = (n+1) % 2;
    
    for(int i = Nx_iono_lower; i < Nx_iono_upper; i++){
        for(int j = Ny_iono_lower; j <= Ny_iono_upper; j++){
            Jex[NEW][i][j] = S[i][j](0,0) * Jex[OLD][i][j] 
                            + S[i][j](0,1) * (Jey[OLD][i+1][j] + Jey[OLD][i][j] + Jey[OLD][i+1][j-1] + Jey[OLD][i][j-1]) / 4.0 
                            + S[i][j](0,2) * ( Jez[OLD][i+1][j] + Jez[OLD][i][j] ) / 2.0 
                            + B[i][j](0,0) * Ex[i][j] 
                            + B[i][j](0,1) * (Ey[i+1][j] + Ey[i][j] + Ey[i+1][j-1] + Ey[i][j-1]) / 4.0 
                            + B[i][j](0,2) * ( Ez[i+1][j] + Ez[i][j] ) / 2.0;
        }
    }
    
    for(int i = Nx_iono_lower; i <= Nx_iono_upper; i++){
        for(int j = Ny_iono_lower; j < Ny_iono_upper; j++){
            Jey[NEW][i][j] = S[i][j](1,0) * (Jex[OLD][i][j+1] + Jex[OLD][i][j] + Jex[OLD][i-1][j+1] + Jex[OLD][i-1][j]) / 4.0 
                            + S[i][j](1,1) * Jey[OLD][i][j]
                            + S[i][j](1,2) * (Jez[OLD][i][j+1] + Jez[OLD][i][j]) / 2.0 
                            + B[i][j](1,0) * (Ex[i][j+1] + Ex[i][j] + Ex[i-1][j+1] + Ex[i-1][j]) / 4.0 
                            + B[i][j](1,1) * Ey[i][j] 
                            + B[i][j](1,2) * (Ez[i][j+1] + Ez[i][j]) / 2.0;
        }
    }

    for(int i = Nx_iono_lower; i <= Nx_iono_upper; i++){
        for(int j = Ny_iono_lower; j <= Ny_iono_upper; j++){
            Jez[NEW][i][j] = S[i][j](2,0) * (Jex[OLD][i][j] + Jex[OLD][i-1][j]) / 2.0 
                            + S[i][j](2,1) * (Jey[OLD][i][j] + Jey[OLD][i][j-1]) / 2.0 
                            + S[i][j](2,2) * Jez[OLD][i][j] 
                            + B[i][j](2,0) * (Ex[i][j] + Ex[i-1][j]) / 2.0 
                            + B[i][j](2,1) * (Ey[i][j] + Ey[i][j-1]) / 2.0 
                            + B[i][j](2,2) * Ez[i][j];
        }
    }

}