#ifndef _utils
#define _utils

extern double energy_to_volume(int, double, double, int, double[][3], int);
extern void calc_corr_function(int, int N, double[], double [][N]);
extern void mean_squared_displacement(int, int N, double [], double [][N], double [][N], double [][N]);
extern void velocity_correlation(int, int N, double [], double [][N], double [][N], double [][N]);
extern void fast_velocity_correlation(int, int N, double [], double [], double [], double [][N], double [][N], double [][N], double);
extern double calc_temp(int, double, double[][3]);
extern double calc_pres(int, double, double[][3], double[][3], double);
extern double calc_temp_scale(double, double, double, double);
extern double calc_pres_scale(double, double, double, double);

#endif
