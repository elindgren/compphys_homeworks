#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#define PI 3.141592653589793238

double trialWavefunction(double alpha, double R1[], double R2[]){
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
    double deltax = x1-x2;
    double deltay = y1-y2;
    double deltaz = z1-z2;


    /* Calculate the energy */
    return -4 + (deltax + deltay + deltaz)/pow((r12 * (1+alpha*r12)),2.0) - 1.0/pow((r12 * (1+alpha*r12)),3.0) - 1.0/(4*pow((r12 * (1+alpha*r12)), 4.0)) + 1.0/r12;
}

double getTheta(double R1[], double R2[]){
    double x1 = R1[0]; double y1 = R1[1]; double z1 = R1[2];
    double x2 = R2[0]; double y2 = R2[1]; double z2 = R2[2];

    /* Calculate r1, r2 and r1*r2 to simplify expression */
    double r1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
    double r2 = sqrt( x2*x2 + y2*y2 + z2*z2 );
    double r1_dot_r2 = x1*x2 + y1*y2 + z1*z2;

    return acos(r1_dot_r2 / (r1*r2));
}

void metropolis(int N, double d, double alpha, double r1[], double r2[], double rho[][3], double theta[], double P_theta[], double E[]){
    /* Performs a MCMC sampling of the configuration space */

    /* Metropolis variables */
    int steps = 0;  // Loop variable
    int acc_steps = 0;
    double P_current;   // Current probability
    double P_proposal;  // Proposed probability
    double P_ratio;
    double *r1_proposal = malloc(3 * sizeof(double));
    double *r2_proposal = malloc(3 * sizeof(double));
    /* RNG */
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
    double u1;  // Random variable
    double u2;  // -||-

    /* Copy positions to proposal arrays */
    for(int i=0; i<3; i++){
        r1_proposal[i] = r1[i];
        r2_proposal[i] = r2[i];
    }

    for(steps=0; steps<N; steps++){
        /* Get current probability */
        P_current = trialWavefunction(alpha, r1, r2)*trialWavefunction(alpha, r1, r2);

        /* Randomize next step - all directions simultanously */
        for(int i=0; i<3; i++){
            u1 = gsl_rng_uniform(q);  // Random value between 0 and 1
            u2 = gsl_rng_uniform(q);  // -||-
            r1_proposal[i] = 2*(u1-0.5) * d;  // Take a randomized, symmetric step with length d
            r2_proposal[i] = 2*(u2-0.5) * d;
        }
        P_proposal = trialWavefunction(alpha, r1_proposal, r2_proposal)*trialWavefunction(alpha, r1_proposal, r2_proposal);

        /* Accept or reject proposal step */
        if(P_proposal / P_current > 1 || P_proposal > gsl_rng_uniform(q)){
            /* Accept the step - overwrite current positions */
            for(int i=0; i<3; i++){
                r1[i] = r1_proposal[i];
                r2[i] = r2_proposal[i];
            }
            acc_steps++;
        }else{
            /* Reject - do nothing */
        }   

        /* Calculate P */
        rho[steps][0] = sqrt( r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2] );  // TODO - calculate correctly
        rho[steps][1] = sqrt( r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2] );  // TODO - calculate correctly
        rho[steps][2] = trialWavefunction(alpha, r1, r2)*trialWavefunction(alpha, r1, r2);

        /* Calculate theta */
        theta[steps] = getTheta(r1, r2);
        P_theta[steps] = 0.5 * sin(theta[steps]);

        /* Calculate E */
        E[steps] = energy(alpha, r1, r2);
    }
    printf("Acceptance ratio over n_steps=%d: %.2f \n", N, (double)acc_steps/steps);

    /* Free variables */
    free(r1_proposal); r1_proposal=NULL; free(r2_proposal); r2_proposal=NULL;
}


void task1(){
    /* Implement the Metropolis algorithm for the helium atom with the trial wave function */

    /****** Initialize variables ******/

    /* Metropolis variables */
    double d = 0.184;  // Displacement parameter - ''step length''
    int n_steps = 10000;  // Number of Metropolis steps
    
    /* System variables */
    double alpha = 0.1;  // Parameter for trial wavefunction
    double *r1 = malloc(3 * sizeof(double));  // Position electron 1
    double *r2 = malloc(3 * sizeof(double));  // Position electron 2
    double *r = malloc(n_steps * sizeof(double));  // Distance from the origin - for rho(r)
    double (*rho)[3] = malloc(sizeof(double[n_steps][3]));  // Sampled probabilities
    double *theta = malloc(n_steps * sizeof(double));  // Sampled values of theta
    double *P_theta = malloc(n_steps * sizeof(double));  // Probability distribution of theta
    double *E = malloc(n_steps * sizeof(double));  // Sampled energies
    /* File handling */
    FILE *f;
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
    metropolis(n_steps, d, alpha, r1, r2, rho, theta, P_theta, E);

    /* Save results */
    f = fopen("task1/rho.dat", "w");
    for(int i=0; i<n_steps; i++){
        for(int j=0; j<3; j++){
            fprintf(f, "%.8f \t", rho[i][j]);   
        }
        fprintf(f, "\n");  
    }
    fclose(f);

    f = fopen("task1/theta.dat", "w");
    for(int i=0; i<n_steps; i++){
        fprintf(f, "%.8f \t %.8f \n", theta[i], P_theta[i]);
    }
    fclose(f);

    f = fopen("task1/energy.dat", "w");
    for(int i=0; i<n_steps; i++){
        fprintf(f, "%.8f \n", E[i]);
    }
    fclose(f);

    /****** Free function variables and fields ******/
    free(r1); r1=NULL; free(r2); r2=NULL; 
    free(rho); rho=NULL; free(theta); theta=NULL;
    gsl_rng_free(q);
}


int main(){
    /* Each task in separate files */
    task1();
}