#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
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
    int NIndepCalc = 10; // Number of independent calculations for each alpha
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
            res[i*NIndepCalc+j][1] = resTup.E;
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
    omp_set_num_threads(4);

    /* Parameters to iterate over */
    int N_steps = 100;  // Number of optimization iterations
    // double alpha;  // Current value of alpha
    double (*res)[2] = malloc(sizeof(double[N_steps][2]));  // Matrix containing results - alpha | E | s_corr | s_block
    /* File IO */
    // FILE *f;

    /***** Iterate over alphas - start each alpha on a separate thread *****/

    /***** Save results *****/
    /***** Free variables *****/
    free(res); res=NULL;
    printf("Finished Task 4 \n");

}