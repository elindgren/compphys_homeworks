#ifndef _utils
#define _utils

extern double calc_temp(int, double, double[][3]);
extern double calc_pres(int, double, double[][3], double[][3], double);
extern double calc_temp_scale(double, double, double, double);
extern double calc_pres_scale(double, double, double, double);

#endif
