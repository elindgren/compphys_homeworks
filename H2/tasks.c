#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include "utils.h"

void runTask1and2(){
    control(0.1, 1);  // Run control for alpha=0.1, and save values according to Task 1 and Task 2
}

void runTask3(){
    /* Run Metropolis for several values of alpha and find optimum value of energy */

    /***** Parameter declarations *****/
    /* OpenMP parameters */
    omp_set_num_threads(8);
    clock_t t; // For calculating CPU time
    double time_taken;
    time_t t_wall;  // For calculating true time taken - ''wall time''


    /* Parameters to iterate over */
    double minAlpha = 0.05;  // Minimum value of alpha
    double maxAlpha = 0.25;  // Maximum value of alpha
    double dAlpha = 0.01;  // Step in alpha
    int N_alphas = (maxAlpha-minAlpha) / dAlpha;  // Number of alphas 
    int NIndepCalc = 100; // Number of independent calculations for each alpha
    double alpha;  // Current value of alpha
    double (*res)[6] = malloc(sizeof(double[N_alphas*NIndepCalc][6]));  // Matrix containing results - alpha | E | s_corr | s_block
    /* File IO */
    FILE *f;

    /***** Iterate over alphas - start each alpha on a separate thread *****/
    printf("Starting Task 3 \n");
    t = clock(); 
    t_wall = time(0);
    #pragma omp parallel for
    for(int i=0; i<N_alphas; i++){
        alpha = minAlpha + i*dAlpha;  
        #pragma omp parallel for
        for(int j=0; j<NIndepCalc; j++){
            // printf("Alpha: %.2f \n", alpha);
            struct resultTuple resTup = control(alpha, 0);  // Declare here to make sure no overwriting of energies between threads
            #pragma omp critical  // This must be executed one thread at a time
            res[i*NIndepCalc+j][0] = alpha;
            res[i*NIndepCalc+j][1] = resTup.meanE;
            res[i*NIndepCalc+j][2] = resTup.varE;
            res[i*NIndepCalc+j][3] = resTup.N;
            res[i*NIndepCalc+j][4] = resTup.sC;
            res[i*NIndepCalc+j][5] = resTup.sB;
        }
    }
    t = clock() - t; 
    t_wall = time(0) - t_wall;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
    printf("\t Finished iteration over alphas - CPU time: %.0f s - Wall time: %ld s\n", time_taken, t_wall);
    /***** Save results *****/
    f = fopen("task3.dat", "w");
        for(int i=0; i<N_alphas*NIndepCalc; i++){
            fprintf(f, "%.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \n", res[i][0], res[i][1], res[i][2], res[i][3], res[i][4], res[i][5]);
        }
    fclose(f);

    /***** Free variables *****/
    free(res); res=NULL;
    printf("Finished Task 3 \n");
}

void runTask4(){
    /* Run Metropolis for several values of alpha and find optimum value of energy using gradient descent */

    /***** Parameter declarations *****/
    /* OpenMP parameters */
    omp_set_num_threads(8);
    clock_t t; // For calculating CPU time
    double time_taken;
    time_t t_wall;  // For calculating true time taken - ''wall time''

    /* Parameters to iterate over */
    double minBeta = 0.5;  // Minimum value of beta
    double maxBeta = 1.0;  // Maximum value of beta
    double dBeta = 0.1;  // Step in beta
    int N_betas = (maxBeta-minBeta) / dBeta;  // Number of beta 
    int N_steps = 1000 + 2; // Number of updates of alpha until convergence - the +2 is the number of header rows
    double (*res)[N_betas] = malloc(sizeof(double[N_steps][N_betas]));  // Matrix containing results - beta[:0] | beta[:1] etc. 
    /* File IO */
    FILE *f;

    /***** Iterate over betas - start each beta on a separate thread *****/
    printf("Starting Task 4 \n");
    t = clock(); 
    t_wall = time(0);
    #pragma omp parallel for
    for(int j=0; j<N_betas; j++){
        /* Declare alpha and beta here for proper scope */
        double alpha = 0.1;  // Starting value of alpha
        double beta = minBeta + j*dBeta;
        res[0][j] = beta;  // Save beta - header
        res[1][j] = alpha;  // Save alpha
        // printf("Beta %.2f \n", beta);
        for(int i=2; i<N_steps; i++){
            struct resultTuple resTup = control(alpha, 0);  // Perform metropolis for current alpha
            alpha -= exp(-beta) * 2*(resTup.meanEGradLnPsi - resTup.meanE*resTup.meanGradLnPsi); // Use returned average energy and value of gradient log to update alpha
            #pragma omp critical  // This must be executed one thread at a time
            res[i][j] = alpha;  // Save alpha
        }
    }
    t = clock() - t; 
    t_wall = time(0) - t_wall;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
    printf("\t Finished iteration over betas - CPU time: %.0f s - Wall time: %ld s\n", time_taken, t_wall);

    /***** Save results *****/
    /* Note that the first row contains all the betas! */
    f = fopen("task4.dat", "w");
        for(int i=0; i<N_steps; i++){
            for(int j=0; j<N_betas; j++){
                fprintf(f, "%.8f \t", res[i][j]);
            }
             fprintf(f, "\n");
        }
    fclose(f);

    /***** Free variables *****/
    free(res); res=NULL;
    printf("Finished Task 4 \n");
}

void runTask5(){
    /* Run Metropolis careful soimulation */

    /***** Parameter declarations *****/
    /* OpenMP parameters */
    omp_set_num_threads(8);
    clock_t t; // For calculating CPU time
    double time_taken;
    time_t t_wall;  // For calculating true time taken - ''wall time''

    /* Parameters to iterate over */
    int N_walkers = 10000;
    double alpha = 0.143; // The optimized value for alpha
    double (*res)[6] = malloc(sizeof(double[N_walkers][6])); 
    int completed_walkers = 0;
    /* File IO */
    FILE *f;

    /***** Iterate over betas - start each beta on a separate thread *****/
    printf("Starting Task 5 \n");
    t = clock(); 
    t_wall = time(0);
    #pragma omp parallel for
    for(int i=0; i<N_walkers; i++){
        if(completed_walkers % 100 == 0){
            printf("Progress: %.2f % \n", (double)completed_walkers/N_walkers*100);
        }
        struct resultTuple resTup = control(alpha, 0);  // Declare here to make sure no overwriting of energies between threads
            #pragma omp critical  // This must be executed one thread at a time
            res[i][0] = alpha;
            res[i][1] = resTup.meanE;
            res[i][2] = resTup.varE;
            res[i][3] = resTup.N;
            res[i][4] = resTup.sC;
            res[i][5] = resTup.sB;
            completed_walkers++;
    }
    t = clock() - t; 
    t_wall = time(0) - t_wall;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
    printf("\t Finished iteration over betas - CPU time: %.0f s - Wall time: %ld s\n", time_taken, t_wall);

    /***** Save results *****/
    /* Note that the first row contains all the betas! */
    f = fopen("task5_test.dat", "w");
        for(int i=0; i<N_walkers; i++){
            fprintf(f, "%.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \n", res[i][0], res[i][1], res[i][2], res[i][3], res[i][4], res[i][5]);
        }
    fclose(f);

    /***** Free variables *****/
    free(res); res=NULL;
    printf("Finished Task 5 \n");
}