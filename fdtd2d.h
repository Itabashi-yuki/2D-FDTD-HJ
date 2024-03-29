#include <cmath>
#include <eigen3/Eigen/Dense>
#include <string>

extern std::string global_dirName;

constexpr double C0 { 3.0e8 };
constexpr double MU0 { 4.0 * M_PI * 1.0e-7 };
constexpr double EPS0 { 1.0 / MU0 / C0 / C0 };
constexpr double CHARGE_e { 1.602e-19 };
constexpr double MASS_e { 9.1e-31 };
constexpr double B0 { 50000e-9 };
// constexpr double B0 { 0.0 };

constexpr double Rx { 100.0e3 };
constexpr double Ry { 100.0e3 };
constexpr double dx { 0.5e3 };
constexpr double dy { 0.5e3 };

constexpr int Nx { int(Rx / dx) };
constexpr int Ny { int(Ry / dy) };

constexpr double Tmax { 0.001 };
// constexpr double Tmax { 0.00085 };
// constexpr double Tmax { 300.0 };

constexpr double Ne_max { 1.0e11 };
constexpr double Omega { std::sqrt(CHARGE_e * CHARGE_e * Ne_max / MASS_e / EPS0) };
constexpr double xi {1.0 / std::sqrt(1 + Omega * Omega / 4.0  / C0 / C0 / ( 1.0 / dx / dx + 1.0 / dy / dy))};
constexpr double dt { xi * 0.9999 / C0 / std::sqrt( 1 / dx / dx + 1 / dy / dy ) };
// constexpr double dt { xi * 1.1 / C0 / std::sqrt( 1 / dx / dx + 1 / dy / dy ) };

// constexpr double Cdt { 0.1 };
// constexpr double dt { Cdt * 0.9 / C0 / std::sqrt( 1 / dx / dx + 1 / dy / dy ) };
constexpr int Nt { int(Tmax / dt) };
// constexpr int Nt { 500 };

constexpr double source_x { 30.0e3 };
constexpr double source_y { 50.0e3 };
// constexpr double source_x { Rx / 2.0 - 20.0e3 };
// constexpr double source_y { Rx / 2.0 };

/* PML parameters */
constexpr double M { 3.7 };
constexpr double R { 1.0e-6 };
constexpr double Z { std::sqrt(MU0 / EPS0) };
constexpr double Y { 1.0 / Z };
constexpr int L { 10 };

/* Plasma parameters */
constexpr double Rx_iono_lower { 50.0e3 };
constexpr double Ry_iono_lower { 40.0e3 };
constexpr double Rx_iono_upper { 70.0e3};
constexpr double Ry_iono_upper { 60.0e3 };
// constexpr double Rx_iono_lower { Rx / 2.0 };
// constexpr double Ry_iono_lower { Ry / 2.0 - 10.0e3 };
// constexpr double Rx_iono_upper { Rx / 2.0 + 20.0e3};
// constexpr double Ry_iono_upper { Ry / 2.0 + 10.0e3 };

constexpr double Nx_iono_lower { int(Rx_iono_lower / dx) };
constexpr double Ny_iono_lower { int(Ry_iono_lower / dy) };
constexpr double Nx_iono_upper { int(Rx_iono_upper / dx) };
constexpr double Ny_iono_upper { int(Ry_iono_upper / dy) };
constexpr double Nx_iono {Nx_iono_upper - Nx_iono_lower + 1};
constexpr double Ny_iono {Ny_iono_upper - Ny_iono_lower + 1};

/* 印加する磁場の方向 */
constexpr double THETA { M_PI / 4.0 }; 
constexpr double PHI { M_PI / 4.0 };

/* 電流源のパラメタ */
// constexpr double f0 { 40.0e3 };
// constexpr double dx_0 { 500 };
// constexpr double dy_0 { 500 };
// constexpr double dt_0 { C0 / std::sqrt( 1.0 / dx_0 / dx_0 + 1.0 / dy_0 / dy_0 ) };
// constexpr double sigma { 70 * dt_0 };
// constexpr double t0 { 4.0 * sigma };
constexpr double sigma { 15 * 3.39e-7 };
constexpr double t0 { 6.0 * sigma };

/* obs parameters */
constexpr double obs_x { 80.0e3 };
constexpr double obs_y { 50.0e3 };
constexpr double obs_t_step { 1.0e-3 };

void update_E(double **Ex, double **Ey, double **Ez, double **Hx, double **Hy, double **Hz,
                 double ***Jex, double ***Jey, double ***Jez, int n);
void update_Ex_PML(double **Ex, double **Hz, double *CEX1, double *CEX2);
void update_Ey_PML(double **Ey, double **Hz, double *CEY1, double *CEY2);
void update_Ez_PML(double **Ez, double **Ezx, double **Ezy, double **Hx, double **Hy,
                    double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2);

void update_H(double **Hx, double **Hy, double **Hz, double **Ex, double **Ey, double **Ez);
void update_Hx_PML(double **Hx, double **Ez, double *CHX1, double *CHX2);
void update_Hy_PML(double **Hy, double **Ez, double *CHY1, double *CHY2);
void update_Hz_PML(double **Hz, double **Hzx, double **Hzy, 
                 double **Ex, double **Ey,double *CHZX1, double *CHZX2,
                 double *CHZY1, double *CHZY2);

void update_J(double ***Jex, double ***Jey, double ***Jez, double **Ex, double **Ey, double **Ez,
                 Eigen::Matrix3d **S, Eigen::Matrix3d **B, int n);
double Jz(double t);

void initialize_PML(double *CEX1, double *CEX2, double *CEY1, double *CEY2,
                    double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2, 
                    double *CHX1, double *CHX2, double *CHY1, double *CHY2,
                    double *CHZX1, double *CHZX2, double *CHZY1, double *CHZY2);
void initialize_Plasma(Eigen::Matrix3d **S, Eigen::Matrix3d **B, double Ne_alt, double nu_alt);

void output_Ez(double **Ex, double **Ey, double **Ez, int n, int n0);
void output_obs(double **Ex, double **Ey, double **Ez, double **Hx, double **Hy, double **Hz,
            double ***Jex, double ***Jey, double ***Jez, int n);
void output_pal(double Ne_alt, double nu_alt);

double cal_nu(double z);
double cal_Ne(double Ne_exp);
int cal_obs_n0();
void judge_Ez_div(double **Ez, int n, int n0, double *Ez_max, double Ne, double nu);

void make_dir();
double *** allocate_3d(int dim1, int dim2, int dim3, double initial_Value);
double ** allocate_2d(int dim1, int dim2, double initial_Value);
double *allocate_1d(int dim1, double initial_Value);
void free_memory3d(double ***array, int dim1, int dim2);
void free_memory2d(double **array, int dim1);