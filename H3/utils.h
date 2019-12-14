#include <gsl/gsl_rng.h>
#ifndef _task1
#define _task1

extern void control(int, int, int, int, double, double, double *, double *);
extern int diffusionMC(int, int, int, double, double *, double *, double, gsl_rng *);
extern double potential(int, int, double *);

#endif