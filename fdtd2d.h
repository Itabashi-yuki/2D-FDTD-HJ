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


void update_E(double **Ez, double **Hx, double **Hy, double *Jz, int n);
void update_E_PML(double **Ez, double **Ezx, double **Ezy, double **Hx, double **Hy, double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2);
void update_H(double **Hx, double **Hy, double **Ez);
void update_H_PML(double **Hx, double **Hy, double **Ez, double *CHX1, double *CHX2, double *CHY1, double *CHY2);
void current_source(double *Jz, int n);
void initialize_PML(double *CEZX1, double *CEZX2, double *CEZY1, double *CEZY2, double *CHX1, double *CHX2, double *CHY1, double *CHY2);
void output_Ez(double **Ez, int n);