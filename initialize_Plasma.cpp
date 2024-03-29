#include "fdtd2d.h"
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <math.h>

constexpr double M2KM {1e-3}; /*[m] to [km]*/

double cal_nu(double z){
    // return 4.303e11 * exp(-0.1622 * z * M2KM);
    double nu_const = std::pow(10,z);
    return nu_const;
    // return 0.0;
}

double omg_p(double Ne){
    return std::sqrt( Ne * CHARGE_e * CHARGE_e / ( MASS_e * EPS0 ) );
}

double cal_Ne(double Ne_exp){
    // double *Alt = new double[Nx];
    // double *Ne = new double[Nx];

    // std::ifstream ifs("Ne.dat");

    // for(int i = Nx_iono_lower; i <= Nx; i++){
    //     ifs >> Alt[i] >> Ne[i];
    // }    
    // constexpr double Ne_input_file_alt_lower = { 60 }; /* 電子密度入力ファイル(Ne.dat)において電子密度が入っている最低高度の値 */

    // /*  入力ファイルの高度と解析領域における位置のズレの修正 */
    // int Ne_alt_int = int(Nx_iono_lower + 2 * (Ne_alt - Ne_input_file_alt_lower));
    // double Ne_const = Ne[Ne_alt_int];
    // double Ne_const = 4.0e5;    
    
    double Ne_const = std::pow(10,Ne_exp);
    // double Ne_const = 0.0;
    return Ne_const;

    // delete [] Alt;
    // delete [] Ne;
}

void initialize_Plasma(Eigen::Matrix3d **S, Eigen::Matrix3d **B, double Ne_exp, double nu_exp){
    double Ne_const = cal_Ne(Ne_exp);
    double nu_const = cal_nu(nu_exp);
    // double nu_const = cal_nu(nu_exp * 1.0e3);

    std::cout << "Ne = " << Ne_const << ", nu = " << nu_const << std::endl;

    for(int i = Nx_iono_lower; i <= Nx_iono_upper; i++){
        for(int j = Ny_iono_lower; j <= Ny_iono_upper; j++){
            double Omg_0 = 1.0 / dt + nu_const / 2.0;
            double Omg_0_prime = 1.0 / dt - nu_const / 2.0;
            double Omg_c = CHARGE_e * B0 / MASS_e;
            double Omg_p = omg_p(Ne_const);

            Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
            Eigen::Matrix3d A, R1, R2;
            Eigen::Vector3d b;

            // const double THETA = M_PI / 2.0;
            // const double PHI = 0.0;

            b << std::sin(THETA) * std::cos(PHI), std::sin(THETA) * std::sin(PHI), std::cos(THETA);

            R2 << 0, b(2), -b(1),
                -b(2), 0, b(0),
                b(1), -b(0), 0;
            
            R1 = Omg_c / 2.0 * R2;

            A = Omg_0 * I + R1;

            S[i][j] = A.inverse() * ( Omg_0_prime * I - R1 );
            B[i][j] = EPS0 * Omg_p * Omg_p * A.inverse();
        }
    }

}