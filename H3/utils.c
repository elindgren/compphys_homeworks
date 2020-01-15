#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define PI 3.141592653589793238


double trial_wave(double x[], double alpha){
    double x1 = x[0];
    double y1 = x[1];
    double z1 = x[2];
    double x2 = x[3];
    double y2 = x[4];
    double z2 = x[5];

    /* Calculate r1, r2 and r1*r2 to simplify expression */
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

    return exp(-2*r1 )*exp(-2*r2) * exp(r12 /(2 * (1+alpha*r12)));
}

double derivative_terms(double x[], int i, double alpha_trial){
    double x1 = x[0];
    double y1 = x[1];
    double z1 = x[2];
    double x2 = x[3];
    double y2 = x[4];
    double z2 = x[5];

    /* Calculate r1, r2 and r1*r2 to simplify expression */
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

    /* Calculate derivative */
    double df1;
    double df2;

    if(i<=2){
        /* Particle 1 - derivatives w.r. xi */
        df1 = (-1/r1) * 2*x[i];
        df2 = (x[i] - x[i+3])/( 2*r12 * pow(1+alpha_trial*r12, 2.0) );
        
        /* From potential! */
        // df1 = -2/pow(r1, 3.0)*x[i];
        // df2 = 1/pow(r12, 3.0) * (x[i] - x[i+3]);
    }else{
        /* Particle 2 - derivative w.r. yi */
        df1 = (-1/r2) * 2*x[i];
        df2 = (x[i] - x[i-3])/( 2*r12 * pow(1+alpha_trial*r12, 2.0) );  // -(x[i]-y[i]) = y[i]-x[i] - xi = x[i-3]
        
        /* From potential! */
        // df1 = -2/pow(r2, 3.0)*x[i];
        // df2 = 1/pow(r12, 3.0) * (x[i] - x[i-3]);
    }
    // printf("df1=%.4f \t df2=%.4f \t r1=%.8f \t r2=%.8f \t r12=%.8f", df1, df2, r1, r2, r12);
    return df1 + df2;
}

double force(double x[], int i){
    /* The trial wavefunction cancels from the gradient, resulting in only 2*derivative factors */
    double alpha_trial = 0.1;
    return 2*derivative_terms(x, i, alpha_trial);
}


double potential(int cs, int ndim, double x[]){
    double V=0;
    if(cs==1){
        for(int k=0; k<ndim; k++){
            V += 0.5*x[k]*x[k];
        }
    }else if(cs==2 || cs==3 || cs==4){
        double x1 = x[0];
        double y1 = x[1];
        double z1 = x[2];
        double x2 = x[3];
        double y2 = x[4];
        double z2 = x[5];

        /* Calculate r1, r2 and r1*r2 to simplify expression */
        double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
        double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
        double r12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
        V = -2/r1 - 2/r2 + 1/r12;
    }
    return V;
}

int diffusionMC(int cs, int ndim, int N, double ET, double x[][ndim], double currentX[], double dtau, gsl_rng *q){
    double V;
    double W;
    double F;
    int m;

    for(int j=0; j<N; j++){
        /* Save current position */
        for(int k=0; k<ndim; k++){
            currentX[k] = x[j][k];  // Pass this to calculate force - wrong particle if not!
        }

        /* Displace walker */
        for(int k=0; k<ndim; k++){
            if(cs==1 || cs==2){
                x[j][k] += sqrt(dtau)*gsl_ran_gaussian(q, 1);
            }else if (cs==3 || cs==4){
                F = force(currentX, k);
                // printf("\t F=%.4f \n", F);
                x[j][k] += sqrt(dtau)*gsl_ran_gaussian(q, 1) + 0.5*dtau*F;
            }
            currentX[k] = x[j][k];
        }

        /* Evaluate weight factor */
        V = potential(cs, ndim, currentX);
        W = exp(dtau * (ET-V));
        // /* Kill or create walkers */  
        m = (int)( W + gsl_rng_uniform(q));
        if(m==0){
            /* Kill this Walker - shift all positions one step to the left */
            for(int l=j+1; l<N; l++){
                for(int k=0; k<ndim; k++){
                    x[l-1][k] = x[l][k];
                }
            }
            N -= 1;
        }else if (m>1){
            /* Create a new walker */
            /* It should spawn m-1 NEW walkers */
            for(int l=N; l<N+m-1; l++){
                for(int k=0; k<ndim; k++){
                    x[l][k] = x[j][k];
                }
            }
            N += m-1; 
        }
        // printf("\t j=%d, m=%d, N=%d, W=%.4f, ET=%.4f, V=%.4f \n", j, m, N, W, ET, V);
        /* Save positions to x_hist */
        
    }
    return N;
}

double evaluateEnergy(int N0, int N, double dtau, double alpha, double ETPhoneHome[], int iteration){
    double avgE = 0;
    for(int l=0; l<iteration; l++){
        avgE += ETPhoneHome[l];
    }
    avgE /= iteration;
    // printf("\t\t AvgE= %.2f logN/N0: %.2f N=%d N0=%d \n", avgE, log((double)N/(double)N0), N, N0);
    
    return avgE - alpha/dtau * log((double)N/(double)N0);
}


void control(int cs, int ndim, int N0, int iters, double dtau, double alpha, double ETPhoneHome[], double Nwalkers[]){
    /****** Initialize variables ******/

    /* Monte Carlo variables */
    int N = N0;
    double ET = -5.0;
    int Nmax = 5*N0;
    double (*x)[ndim] = malloc(sizeof(double[Nmax][ndim]));
    double *currentX = malloc(ndim*sizeof(double));

    int bins = 500;
    double xmax = 5;
    int (*x_hist)[bins] = malloc(sizeof(double[iters][bins]));

    /* RNG */
    const gsl_rng_type *T;
    gsl_rng *q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q, time(NULL));

    /* Randomize x starting values */
    for(int l=0; l<N0; l++){
        for(int k=0; k<ndim; k++){
            x[l][k] = gsl_ran_gaussian(q, 1);
        }
    }

    /* Save initial ET and N */
    Nwalkers[0] = N;
    ETPhoneHome[0] = ET;
    
    for(int i=1; i<iters+1; i++){
        // printf("Iteration: %d \n", i);
        N = diffusionMC(cs, ndim, N, ET, x, currentX, dtau, q);
        // printf("\t N=%d \n", N);
        ET = evaluateEnergy(N0, N, dtau, alpha, ETPhoneHome, i);
        // printf("\t ET=%.2f \n", ET);

        double pos;
        int idx;
        if(cs==1){
            for(int j=0; j<N; j++){
                pos = x[j][0];
                idx = (int)(bins/2 * (pos/xmax + 1));
                if(idx < 0){
                    idx = 0;
                }else if(idx > bins-1){
                    idx = bins-1;
                }
                x_hist[i][idx] += 1; 
            }
        }
        /* Save results */
        Nwalkers[i] = N;
        ETPhoneHome[i] = ET;
    }

    if(cs==1){
        /* Write x_hist to file */
        FILE *f;

        /***** Write to file *****/
        double tau;
        f = fopen("task1_hist.dat", "w");
            for(int i=0; i<iters; i++){
                tau = i*dtau;
                fprintf(f, "%.8f \t", tau);
                for(int j=0; j<bins; j++){
                     fprintf(f, "%d \t", x_hist[i][j]);
                }
                fprintf(f, "\n");
            }
        fclose(f);
    }

    /****** Free variables ******/
    free(x); x=NULL;
    free(currentX); currentX=NULL;
    free(x_hist); x_hist=NULL;
}


