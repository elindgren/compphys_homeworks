#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "utils.h"
#define PI 3.141592653589793238

double trialWavefunction(double alpha, double R1[], double R2[])
{
    double x1 = R1[0];
    double y1 = R1[1];
    double z1 = R1[2];
    double x2 = R2[0];
    double y2 = R2[1];
    double z2 = R2[2];

    /* Calculate r1, r2 and r12 to simplify expression */
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

    /* Calculate and return wave function */
    return exp(-2 * r1) * exp(-2 * r2) * exp(r12 / (2 * (1 + alpha * r12)));
}

double getProb(double alpha, double R1[], double R2[])
{
    return trialWavefunction(alpha, R1, R2) * trialWavefunction(alpha, R1, R2);
}

double getEnergy(double alpha, double R1[], double R2[])
{
    double x1 = R1[0];
    double y1 = R1[1];
    double z1 = R1[2];
    double x2 = R2[0];
    double y2 = R2[1];
    double z2 = R2[2];

    /* Calculate r1, r2 and r12 to simplify expression */
    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
    double r1Dotr2 = x1*x2 + y1*y2 + z1*z2;
    
    /* Define energy terms */
    double E0 = -4.0;
    double E1 = (r1 + r2 - r1Dotr2*(1.0/r1 + 1.0/r2)) / (r12 * pow(1+alpha*r12, 2.0));
    double E2 = -1.0 / (r12 * pow(1+alpha*r12, 3.0));
    double E3 = -1.0 / (4 * pow(1+alpha*r12, 4.0));
    double E4 = 1.0 / r12;

    /* Calculate the energy */
    // -4.0 + (r1 + r2 - r1Dotr2 * (1 / r1 + 1 / r2)) / (r12 * pow(1 + alpha * r12, 2.0)) - 1.0 / (r12 * pow(1 + alpha * r12, 3.0)) - 1.0 / (4 * pow(1 + alpha * r12, 4.0)) + 1.0 / r12;
    return E0 + E1 + E2 + E3 + E4;
}

double getTheta(double R1[], double R2[])
{
    double x1 = R1[0];
    double y1 = R1[1];
    double z1 = R1[2];
    double x2 = R2[0];
    double y2 = R2[1];
    double z2 = R2[2];

    /* Calculate r1, r2 and r1*r2 to simplify expression */
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r1_dot_r2 = x1 * x2 + y1 * y2 + z1 * z2;

    return acos(r1_dot_r2 / (r1 * r2));
}

double getGradLnPsi(double alpha, double R1[], double R2[]){
    double x1 = R1[0];
    double y1 = R1[1];
    double z1 = R1[2];
    double x2 = R2[0];
    double y2 = R2[1];
    double z2 = R2[2];

    /* Calculate r12 to simplify expression */
    double r12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

    return -pow(r12,2.0) / (2 * pow(1+alpha*r12, 2.0) );
}

/* Calculate the statistical inefficiency for the energy based on the correlation function method truncated at kMax */
int statInCorrMethod(double Phi[], double E[], int N, int kMax)
{
    int i;
    int k;
    int s = 0;
    double sum;
    double meanE = 0;   // Mean energy value
    double meanESq = 0; // Mean squared energy value

    /* Pre-calculate the mean and squared mean needed for the correlation functions */
    for (i = 0; i < N; i++)
    {
        meanE += E[i];
        meanESq += pow(E[i], 2.0);
    }
    meanE /= N;
    meanESq /= N;
    //printf("Data statistics for testing: %e %e %e \n", meanE, meanESq, meanESq - pow(meanE, 2.0));

    /* Calculate correlation function up to kMax */
    // printf("Correlation function: \n");
    for (k = 0; k < kMax; k++)
    {
        sum = 0.0;
        for (i = 0; i < N - k; i++)
        {
            sum += E[i] * E[i + k] - pow(meanE, 2.0);
        }
        Phi[k] = sum / (N - k);
        Phi[k] /= meanESq - pow(meanE, 2.0);
        // printf("%d %.6f \n", k, Phi[k]);
    }

    /* Calculate the statistical inefficiency s by finding the index corresponding to the first index for Phi <= exp(-2) */
    for (k = 0; k < kMax; k++)
    {
        if (Phi[k] <= exp(-2.0))
        {
            s = k;
            break;
        }
    }
    // printf("Statistical inefficiency: %d \n", s);
    return s;
}

/* Calculate the statistical inefficiency for the energy based on the block averaging method truncated at BMax */
double statInBlockMethod(double S[], double E[], int N, int BMax)
{
    int i;
    int j;
    int B;
    double s = 0.0;
    double sum;
    double meanE = 0;   // Mean energy value
    double meanESq = 0; // Mean squared energy value
    double meanF;       // Block average
    double meanFSq;     // Squared block average

    /* Pre-calculate the mean and squared mean needed for the correlation functions */
    for (i = 0; i < N; i++)
    {
        meanE += E[i];
        meanESq += pow(E[i], 2.0);
    }
    meanE /= N;
    meanESq /= N;
    // printf("Data statistics for testing: %e %e %e \n", meanE, meanESq, meanESq - pow(meanE, 2.0));

    /* Calculate the lower bound of s as a function of B up to BMax */
    // printf("Lower bound statistical inefficiency: \n");
    for (B = 1; B <= BMax; B++)
    {
        meanF = 0.0;
        meanFSq = 0.0;
        for (j = 0; j < N / B; j++) // Note integer division
        {
            sum = 0.0;
            for (i = 0; i < B; i++)
            {
                sum += E[i + j * B];
            }
            sum /= B;
            meanF += sum;
            meanFSq += pow(sum, 2.0);
        }
        meanF /= N / B; // Integer division with number of block partitions
        meanFSq /= N / B;
        S[B - 1] = B * (meanFSq - pow(meanF, 2.0)) / (meanESq - pow(meanE, 2.0));
    }

    /* Calculate the statistical inefficiency s as the mean of the last 10 % of the evaluated lower bounds */
    for (i = (int)(0.9 * BMax); i < BMax; i++)
    {
        // printf("%d %.6f \n", i, S[i]);
        s += S[i];
    }
    s /= (int)(0.1 * BMax);

    // printf("Statistical inefficiency: %f \n", s);
    return s;
}

void metropolis(int N, double d, double alpha, double r1[], double r2[], double rho[][3], double theta[], double P_theta[], double E[], double gradLnPsi[], int task1and2)
{
    /* Performs a MCMC sampling of the configuration space */

    /* Metropolis variables */
    int steps = 0; // Loop variable
    int acc_steps = 0;
    double P_current;  // Current probability
    double P_proposal; // Proposed probability
    double P_accept;
    double *r1_proposal = malloc(3 * sizeof(double));
    double *r2_proposal = malloc(3 * sizeof(double));
    /* RNG */
    const gsl_rng_type *T;
    gsl_rng *q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q, time(NULL));
    double u; // Random variable

    for (steps = 0; steps < N; steps++)
    {
        /* Get current probability */
        P_current = getProb(alpha, r1, r2);

        /* Copy positions to proposal arrays */
        for (int i = 0; i < 3; i++)
        {
            r1_proposal[i] = r1[i];
            r2_proposal[i] = r2[i];
        }

        /* Randomize next step - all directions simultanously */
        for (int i = 0; i < 3; i++)
        {
            r1_proposal[i] += 2 * (gsl_rng_uniform(q) - 0.5) * d; // Take a randomized, symmetric step with length d
            r2_proposal[i] += 2 * (gsl_rng_uniform(q) - 0.5) * d;
        }
        P_proposal = getProb(alpha, r1_proposal, r2_proposal);

        /* Accept or reject proposal step */
        P_accept = P_proposal / P_current;
        u = gsl_rng_uniform(q);
        // printf("u: %.4f \n", u);
        // printf("P_current: %.4f \n", P_current);
        // printf("P_accept: %.4f \n", P_accept);
        if (P_accept > u)
        {
            /* Accept the step - overwrite current positions */
            // printf("Accepted\n\n");
            for (int i = 0; i < 3; i++)
            {
                r1[i] = r1_proposal[i];
                r2[i] = r2_proposal[i];
            }
            acc_steps++;
        }
        else
        {
            /* Reject - do nothing */
            // printf("Rejected\n\n");
        }
        
        if (task1and2)
        {
            /* Calculate P */
            rho[steps][0] = sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
            rho[steps][1] = sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]);
            rho[steps][2] = getProb(alpha, r1, r2);

            /* Calculate theta */
            theta[steps] = getTheta(r1, r2);
            P_theta[steps] = 0.5 * sin(theta[steps]);
        }
        /* Calculate E */
        E[steps] = getEnergy(alpha, r1, r2);
        /* Calculate grad{ln{Psi_T}} */
        gradLnPsi[steps] = getGradLnPsi(alpha, r1, r2);
     
    }

    if (task1and2)
    {
        printf("Acceptance ratio over N=%d steps: %.4f \n", N, (double)acc_steps / steps);
    }
    /* Free variables */
    free(r1_proposal);
    r1_proposal = NULL;
    free(r2_proposal);
    r2_proposal = NULL;
}

struct resultTuple control(double alpha, int task1and2)
{
    /* Implement the Metropolis algorithm for the helium atom with the trial wave function */

    /****** Initialize variables ******/

    /* Metropolis variables */
    // double d = 0.68;  // Displacement parameter - ''step length''
    double d = 0.485;
    int kMax = 200;       // Maximum number of correlation function evaluations
    int BMax = 2000;      // Maximum number of block length. Must be less than N_tot
    int N_tot = 75000;   // Number of Metropolis steps
    int N_eq = 500;       // Number of equilibration steps - Set this to 0 for task1
    int N = N_tot - N_eq; // Number of production steps (return variable)

    /* System variables */
    double meanE = 0;                                   // Mean energy (return variable)
    double varE = 0;                                    // Variance in energy (return variable)
    int sC = 0;                                         // Statistical ineffieicy correlation method (return variable)
    double sB = 0;                                      // Statistical ineffieicy block averaging method (return variable)
    double meanGradLnPsi = 0;                           // Mean value of gradient of logarithm of wavefunction
    double meanEGradLnPsi = 0;                          // Mean value of Energy and gradient of logarithm of wavefunction
    struct resultTuple resTup = {0, 0, 0, 0, 0};        // Create and initialize results tuple
    double *r1 = malloc(3 * sizeof(double));            // Position electron 1
    double *r2 = malloc(3 * sizeof(double));            // Position electron 2
    double(*rho)[3] = malloc(sizeof(double[N_tot][3])); // Sampled probabilities
    double *theta = malloc(N_tot * sizeof(double));     // Sampled values of theta
    double *P_theta = malloc(N_tot * sizeof(double));   // Probability distribution of theta
    double *E = malloc(N_tot * sizeof(double));         // Sampled energies
    double *Phi = malloc(kMax * sizeof(double));        // Correlation function
    double *S = malloc(BMax * sizeof(double));          // Lower bound of statistical inefficiency s
    double *gradLnPsi = malloc(N_tot * sizeof(double));
    /* File handling */
    FILE *f;
    /* RNG */
    const gsl_rng_type *T;
    gsl_rng *q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q, time(NULL));

    /****** Initialize positions for both electrons ******/
    for (int i = 0; i < 3; i++)
    {
        r1[i] = 2 * (gsl_rng_uniform(q) - 0.5) * 10; // TODO arbitrary as of now
        r2[i] = 2 * (gsl_rng_uniform(q) - 0.5) * 10; // TODO arbitrary as of now
    }

    /****** Metropolis ******/
    metropolis(N_tot, d, alpha, r1, r2, rho, theta, P_theta, E, gradLnPsi, task1and2);

    /****** Task 2 - Statistical inneficiency ******/
    /* Calculate the statistical inneficiency in the sampled energies */
    sC = statInCorrMethod(Phi, E, N_tot, kMax);
    sB = statInBlockMethod(S, E, N_tot, BMax);

    /****** Save results ******/
    if (task1and2)
    {
        /* Save to file if Task 1 or Task 2*/
        f = fopen("task1and2/rho.dat", "w");
        for (int i = 0; i < N_tot; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                fprintf(f, "%.8f \t", rho[i][j]);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        f = fopen("task1and2/theta.dat", "w");
        for (int i = 0; i < N_tot; i++)
        {
            fprintf(f, "%.8f \t %.8f \n", theta[i], P_theta[i]);
        }
        fclose(f);

        f = fopen("task1and2/energy.dat", "w");
        for (int i = 0; i < N_tot; i++)
        {
            fprintf(f, "%.8f \n", E[i]);
        }
        fclose(f);
    }
    else
    {
        /* Return mean value of E and the statistical inefiency */
        for (int i = N_eq; i < N_tot; i++)
        {
            meanE += E[i];
            varE += E[i] * E[i];
            meanGradLnPsi += gradLnPsi[i];
            meanEGradLnPsi += E[i]*gradLnPsi[i];
        }
        meanE /= N;
        varE /= N;             // Mean squared value
        varE -= meanE * meanE; // subtract mean squared
        meanGradLnPsi /= N;
        meanEGradLnPsi /= N;
        // sB = 0;
        // sC = 0;
    }

    /****** Free function variables and fields ******/
    free(r1);
    r1 = NULL;
    free(r2);
    r2 = NULL;
    free(rho);
    rho = NULL;
    free(theta);
    theta = NULL;
    free(P_theta); P_theta = NULL;
    free(E); E = NULL;
    free(Phi); Phi = NULL;
    free(S); S = NULL;
    free(gradLnPsi); gradLnPsi = NULL;
    gsl_rng_free(q);

    /* Return E and s */
    resTup.meanE = meanE;
    resTup.varE = varE;
    resTup.N = N;
    resTup.sB = sB;
    resTup.sC = sC;
    resTup.meanGradLnPsi = meanGradLnPsi;
    resTup.meanEGradLnPsi = meanEGradLnPsi;
    return resTup;
}