#include <cmath>

constexpr double C0 { 3.0e8 };
constexpr double MU0 { 4.0 * M_PI * 1e-7 };
constexpr double EPS0 { 1.0 / MU0 / C0 / C0 };

constexpr double Rx { 100.0e3 };
constexpr double Ry { 100.0e3 };
constexpr double dx { 1.0e3 };
constexpr double dy { 1.0e3 };

constexpr int Nx { int(Rx / dx) };
constexpr int Ny { int(Ry / dy) };

constexpr double Tmax { 0.00085 };
constexpr double dt { 0.9 / C0 / std::sqrt( 1 / dx / dx + 1 / dy / dy ) };
constexpr int Nt { int(Tmax / dt) };

/* PML parameters */
constexpr double M { 3.7 };
constexpr double R { 1.0e-6 };
constexpr double Z { std::sqrt(MU0 / EPS0) };
constexpr double Y { 1.0 / Z };
constexpr int L { 10 };

// struct PML_idx{
//     int i1, j1, i2, j2;
// };

// class PML_domain_E{
// public:
//     PML_idx idx[4];

//     PML_domain_E(){
//         idx[0] = {1, 1, L, Ny - 1};
//         idx[1] = {Nx - L, 1, Nx - 1, Ny - 1};
//         idx[2] = {L + 1, 1, Nx - L - 1, L};
//         idx[3] = {L + 1, Ny - L, Nx - L - 1, Ny - 1};
//     }      
// };

// class PML_domain_Hx{
// public:
//     PML_idx idx[4];

//     PML_domain_Hx(){
//         idx[0] = {1, 0, L, Ny - 1};
//         idx[1] = {Nx - L, 0, Nx - 1, Ny - 1};
//         idx[2] = {L + 1, 0, Nx - L - 1, L - 1};
//         idx[3] = {L + 1, Ny - L, Nx - L - 1, Ny - 1};
//     }      
// };

// class PML_domain_Hy{
// public:
//     PML_idx idx[4];

//     PML_domain_Hy(){
//         idx[0] = {0, 1, L - 1, Ny - 1};
//         idx[1] = {Nx - L, 1, Nx - 1, Ny - 1};
//         idx[2] = {L, 1, Nx - L - 1, L};
//         idx[3] = {L, Ny - L, Nx - L - 1, Ny - 1};
//     }
// };

void update_E(double **Ez, double **Hx, double **Hy, double *Jz, int n);
void update_E_PML(double **Ez, double **Ezx, double **Ezy, double **Hx, double **Hy, double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2);
void update_H(double **Hx, double **Hy, double **Ez);
void update_H_PML(double **Hx, double **Hy, double **Ez, double *CHX1, double *CHX2, double *CHY1, double *CHY2);
void current_source(double *Jz, int n);
void initialize_PML(double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2, double *CHX1, double *CHX2, double *CHY1, double *CHY2);
void output_Ez(double **Ez, int n);