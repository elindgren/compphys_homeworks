#ifndef _task1
#define _task1

struct resultTuple
{
    double E;  // Average energy
    double sC; // s using correlation function
    double sB; // s using block average
};

extern double getProb(double, double[], double[]);
extern double getEnergy(double, double[], double[]);
extern double getTheta(double[], double[]);
extern int statInCorrMethod(double[], double[], int, int);
extern void statInBlockMethod(double[], double[], int, int);
extern void metropolis(int, double, double, double[], double[], double[][3], double[], double[], double[], int);
extern struct resultTuple control(double, int);

#endif