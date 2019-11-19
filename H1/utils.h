#ifndef _utils
#define _utils

extern double energy_to_volume(int, double, double, int, double[][3], int);
extern void calc_corr_function(int N, double *, double[][N], int, int);
extern double calc_temp(int, double, double[][3]);
extern double calc_pres(int, double, double[][3], double[][3], double);
extern double calc_temp_scale(double, double, double, double);
extern double calc_pres_scale(double, double, double, double);

#endif
