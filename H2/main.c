#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#define PI 3.141592653589793238

double trial_wavefunction(double alpha, double R1[], double R2[]){
    double x1 = R1[0]; double y1 = R1[1]; double z1 = R1[2];
    double x2 = R2[0]; double y2 = R2[1]; double z2 = R2[2];

    /* Calculate r1, r2 and r12 to simplify expression */
    double r1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
    double r2 = sqrt( x2*x2 + y2*y2 + z2*z2 );
    double r12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

    /* Calculate and return wave function */
    return exp(-2 * r1) * exp(-2 * r2) * exp( r12 / (2 * (1 + alpha*r12)));
}

double energy(double alpha, double R1[], double R2[]){
    double x1 = R1[0]; double y1 = R1[1]; double z1 = R1[2];
    double x2 = R2[0]; double y2 = R2[1]; double z2 = R2[2];

    /* Calculate r1, r2 and r12 to simplify expression */
    double r1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
    double r2 = sqrt( x2*x2 + y2*y2 + z2*z2 );
    double r12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

    /* Calculate the energy */

    return 0;

}


void task1(){
    /* Implement the Metropolis algorithm for the helium atom with the trial wave function */

    /****** Initialize variables ******/

    /* Metropolis variables */
    double d = 1;  // Displacement parameter - ''step length''
    int n_steps = 10;  // Number of Metropolis steps
    double P_current;   // Current probability
    double P_proposal;  // Proposed probability
    /* System variables */
    double alpha = 0.1;  // Parameter for trial wavefunction
    double *r1 = malloc(3 * sizeof(double));  // Position electron 1
    double *r2 = malloc(3 * sizeof(double));  // Position electron 2
    double *rho = malloc(n_steps * sizeof(double));  // Sampled probabilities
    double *theta = malloc(n_steps * sizeof(double));  // Sampled values of theta
    /* RNG */
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
    double u;  // Random variable

    /****** Initialize positions for both electrons ******/
    for(int i = 0; i<3; i++){
        r1[i] = gsl_rng_uniform(q) * 1; // TODO arbitrary as of now
        r2[i] = gsl_rng_uniform(q) * 1; // TODO arbitrary as of now
    }

    /****** Metropolis ******/
    

    /****** Free function variables and fields ******/
    free(r1); r1=NULL; free(r2); r2=NULL; 
    free(rho); rho=NULL; free(theta); theta=NULL;
    gsl_rng_free(q);
}


int main(){
    /* Each task in separate files */
    task1();
}