#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include "utils.h"

void runTask1(){
    /***** Parameter declarations *****/
    int iters = 20000;  // Monte Carlo iterations
    int cs = 1;  // cs = case för att c suger röv
    int ndim = 1;
    int N0 = 1000;
    double dtau = 0.05;
    double alpha = 0.1;
    /* Data structures */
    double *Nwalkers = malloc(iters * sizeof(double));
    double *ETPhoneHome = malloc(iters * sizeof(double));
    /* File IO */
    FILE *f;

    /***** Run diffusion Monte Carlo *****/
    control(cs, ndim, N0, iters, dtau, alpha, ETPhoneHome, Nwalkers);  // Run control for alpha=0.1, and save values according to Task 1 and Task 2

    /***** Write to file *****/
    double tau;
    f = fopen("task1.dat", "w");
        for(int i=0; i<iters; i++){
            tau = i*dtau;
            fprintf(f, "%.8f \t %.8f \t %.8f \n", tau, Nwalkers[i], ETPhoneHome[i]);
        }
    fclose(f);

    /***** Free variables to get that sweet memory back *****/
    free(ETPhoneHome); ETPhoneHome=NULL;
    free(Nwalkers); Nwalkers=NULL;
}

void runTask2(){
    /***** Parameter declarations *****/
    int iters = 5000;  // Monte Carlo iterations
    int cs = 2;  // cs = case för att c suger röv
    int ndim = 6;
    int N0 = 500;
    double dtau = 0.01;
    double alpha = 0.25;
    /* Data structures */
    double *Nwalkers = malloc(iters * sizeof(double));
    double *ETPhoneHome = malloc(iters * sizeof(double));
    /* File IO */
    FILE *f;

    /***** Run diffusion Monte Carlo *****/
    control(cs, ndim, N0, iters, dtau, alpha, ETPhoneHome, Nwalkers);  // Run control for alpha=0.1, and save values according to Task 1 and Task 2

    /***** Write to file *****/
    double tau;
    f = fopen("task2.dat", "w");
        for(int i=0; i<iters; i++){
            tau = i*dtau;
            fprintf(f, "%.8f \t %.8f \t %.8f \n", tau, Nwalkers[i], ETPhoneHome[i]);
        }
    fclose(f);

    /***** Free variables to get that sweet memory back *****/
    free(ETPhoneHome); ETPhoneHome=NULL;
    free(Nwalkers); Nwalkers=NULL;
}

void runTask3(){
    /***** Parameter declarations *****/
    int iters = 5000;  // Monte Carlo iterations
    int cs = 3;  // cs = case för att c suger röv
    int ndim = 6;
    int N0 = 500;
    double dtau = 0.005;
    double alpha = 0.05;
    /* Data structures */
    double *Nwalkers = malloc(iters * sizeof(double));
    double *ETPhoneHome = malloc(iters * sizeof(double));
    /* File IO */
    FILE *f;

    /***** Run diffusion Monte Carlo *****/
    control(cs, ndim, N0, iters, dtau, alpha, ETPhoneHome, Nwalkers);  // Run control for alpha=0.1, and save values according to Task 1 and Task 2

    /***** Write to file *****/
    double tau;
    f = fopen("task3.dat", "w");
        for(int i=0; i<iters; i++){
            tau = i*dtau;
            fprintf(f, "%.8f \t %.8f \t %.8f \n", tau, Nwalkers[i], ETPhoneHome[i]);
        }
    fclose(f);

    /***** Free variables to get that sweet memory back *****/
    free(ETPhoneHome); ETPhoneHome=NULL;
    free(Nwalkers); Nwalkers=NULL;
}
