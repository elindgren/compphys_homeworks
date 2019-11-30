/*
 MD_main2.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "initfcc.h"
#include "alpotential.h"
#include "utils.h"
#include "fft_func.h"
#define N 256              // Number of atoms
#define eq_timesteps 20000 // Number of equilibration timesteps
#define prod_timesteps 10000   // Number of production

/*
 * Encapsulates the velocity Verlet algorithm
 * Lattice constant a, number of dimensions, number of unit cells in all directions and mass of Al
 * Modified for task 5 to calculate the mean squared displacement as a function of time
 */

double calc_ek(double v[][3], double m_al)
{
    /* Calculate the kinetic energy - T = 0.5*m*|v|^2 */
    double ek = 0;
    double abs_v_sq = 0;
    for (int j = 0; j < N; j += 1)
    {
        abs_v_sq = 0;
        for (int k = 0; k < 3; k += 1)
        {
            abs_v_sq += v[j][k] * v[j][k]; // ek = 0.5*mv^2
        }
        ek += 0.5 * m_al * abs_v_sq;
    }
    return ek;
}

void control(double x[][3], double v[][3], double a[][3], double F[][3], double *a_lat, int ndim, int Nc, double dt, double m_al, int equilibrate, double Teq, double Peq, char label[])
{
    /* Declarations */
    int i, j, k; // loop variables

    int timesteps;
    if (equilibrate)
    {
        timesteps = eq_timesteps;
    }
    else
    {
        timesteps = prod_timesteps;
    }

    double(*X)[N] = malloc(sizeof(double[timesteps + 1][N])); // Position in x-dimension
    double(*Y)[N] = malloc(sizeof(double[timesteps + 1][N]));
    double(*Z)[N] = malloc(sizeof(double[timesteps + 1][N]));

    double(*Vx)[N] = malloc(sizeof(double[timesteps + 1][N])); // Velocity in y-dimension
    double(*Vy)[N] = malloc(sizeof(double[timesteps + 1][N]));
    double(*Vz)[N] = malloc(sizeof(double[timesteps + 1][N]));

    double *Ep = malloc(sizeof(double[timesteps + 1])); // Potential energy
    double *Ek = malloc(sizeof(double[timesteps + 1])); // Kinetic energy
    double *Et = malloc(sizeof(double[timesteps + 1])); // Total energy

    /* File IO */
    FILE *f;
    double t;
    char filename[50] = "";
    char dir[100] = "datafiles/";

    /* Task 3 */
    double T, P;                      // Pressure and temperature
    double alpha_t, alpha_p;          // Scaling parameters for equilibration
    double timedecay_temp = 250 * dt; // Timedecay for temperature scaling
    double timedecay_pres = 500 * dt; // Timedecay for pressure scaling
    double *Temperatures = malloc(sizeof(double[timesteps + 1]));
    double *Pressures = malloc(sizeof(double[timesteps + 1]));
    double *Lat_params = malloc(sizeof(double[timesteps + 1]));

    int p = 0; // Particle to save position for

    /* Task 5 */
    double *MSD = malloc((timesteps + 1) * sizeof(double)); // Mean squared displacement as measured from the start point

    /* Task 6 */
    double *Phi = malloc(sizeof(double[timesteps + 1])); // Time velocity correlation function - not defined for last corr_offset values.

    /* Task 7 */
    double *freq = malloc(2 * (timesteps + 1) * sizeof(double));
    double *Powerspectrum = malloc(2 * (timesteps + 1) * sizeof(double));
    double *fast_Phi = malloc(2 * (timesteps + 1) * sizeof(double));

    /* Calculate energies for initial conditions */
    Ep[0] = get_energy_AL(x, Nc * *a_lat, N); // Supercell length is Nc*a_lat
    Ek[0] = calc_ek(v, m_al);
    Et[0] = Ep[0] + Ek[0]; // Total energy is sum of kinetic and potential

    /* Save initial positions and velocities for particles */
    for (j = 0; j < N; j += 1)
    {
        X[0][j] = x[j][0];
        Y[0][j] = x[j][1];
        Z[0][j] = x[j][2];

        Vx[0][j] = v[j][0];
        Vy[0][j] = v[j][1];
        Vz[0][j] = v[j][2];
    }

    /* The velocity Verlet algorithm */
    for (i = 1; i < timesteps + 1; i += 1)
    {
        // velocity_verlet_step(x, v, a, F, m_al, dt, Nc, &a_lat);
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < 3; k += 1)
            {
                v[j][k] += dt / 2 * a[j][k];
            }
        }

        /* q(t+dt) */
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < 3; k += 1)
            {
                x[j][k] += dt * v[j][k];
            }
        }

        /* a(t+dt) */
        get_forces_AL(F, x, Nc * *a_lat, N); // Supercell length is Nc*a_lat
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < 3; k += 1)
            {
                a[j][k] = F[j][k] / m_al;
                // printf("a[0]: %.2f \t a[1]: %.2f \t  a[2]: %.2f \n", a[j][0], a[j][1], a[j][2]);
                // a[j][k] = 0.0;
            }
        }

        /* v(t+dt) */
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < 3; k += 1)
            {
                v[j][k] += dt / 2 * a[j][k];
            }
        }

        /* Task 2 - Calculate and save energies for this iteration */
        Ep[i] = get_energy_AL(x, Nc * *a_lat, N); // Supercell length is Nc*a_lat
        Ek[i] = calc_ek(v, m_al);
        Et[i] = Ep[i] + Ek[i]; // Total energy is sum of kinetic and potential

        /* Task 3 - Save positions and velocities */
        for (j = 0; j < N; j++)
        {
            X[i][j] = x[j][0];
            Y[i][j] = x[j][1];
            Z[i][j] = x[j][2];

            Vx[i][j] = v[j][0];
            Vy[i][j] = v[j][1];
            Vz[i][j] = v[j][2];
        }

        /* Set scaling parameters */
        T = calc_temp(N, m_al, v);
        Temperatures[i] = T;
        P = calc_pres(N, m_al, v, x, Nc * *a_lat);
        Pressures[i] = P;
        if (equilibrate)
        {
            alpha_t = calc_temp_scale(dt, timedecay_temp, Teq, T);
            alpha_p = calc_pres_scale(dt, timedecay_pres, Peq, P);
        }
        else
        {
            alpha_t = 1;
            alpha_p = 1;
        }

        /* Rescale for equilibration */
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < ndim; k += 1)
            {
                v[j][k] = v[j][k] * sqrt(alpha_t);           // Rescales velocity - changes temperature
                x[j][k] = x[j][k] * pow(alpha_p, 1.0 / 3.0); // Rescales positions
            }
        }
        double a_temp = *a_lat * pow(alpha_p, 1.0 / 3.0); // Rescales volume - changes pressure
        *a_lat = a_temp;
        Lat_params[i] = *a_lat;
    }

    /* Calculate velocity correlation function */
    /* Un-comment to include calculations
    if (!equilibrate)
    {   
        omp_set_num_threads(8);
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                printf("Calculating MSD function \n");
                mean_squared_displacement(timesteps+1, N, MSD, X, Y, Z);  // Task 5 
                printf("Finished MSD function \n");
            }
            #pragma omp section
            {
                printf("Calculating velocity correlation function \n");
                velocity_correlation(timesteps+1, N, Phi, Vx, Vy, Vz); // Task 6
                printf("Finished velocity correlation function \n");
            }
            #pragma omp section
            {
                printf("Calculating power spectrum function \n");
                fast_velocity_correlation(timesteps+1, N, fast_Phi, Powerspectrum, freq, Vx, Vy, Vz, dt);  // Task 7 
                printf("Finished power spectrum function \n");
            }
        }
    }
    */

    strcat(dir, label);
    strcat(filename, dir);
    if (equilibrate)
    {
        /* Equilibration run */
        printf("Saving equilibration \n");

        filename[0] = '\0'; // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/equilibration.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.8f \t %.8f \t %.8f \t %.8f \n", t, Temperatures[i], Pressures[i], Lat_params[i]);
        }
        fclose(f);
    }
    else
    {
        /* Production run */
        printf("Saving production \n");

        // /* Save temperature, pressure and latice constant */
        filename[0] = '\0'; // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/temp_pres_lat.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.8f \t %.8f \t %.8f \t %.8f \n", t, Temperatures[i], Pressures[i], Lat_params[i]);
        }
        fclose(f);

        // /* Save energies to file */
        filename[0] = '\0'; // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/vv_energies.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.4f \t %.4f \t %.4f \t %.4f \n", t, Ep[i], Ek[i], Et[i]);
        }
        fclose(f);

        /* Save position of particle p and MSD to file */
        filename[0] = '\0'; // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/MSD.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", t, X[i][p], Y[i][p], Z[i][p], MSD[i]);
        }
        fclose(f);

        /* Save velocity correlation function to file */
        filename[0] = '\0'; // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/vel_corr.dat");
        f = fopen(filename, "w");
        t = 0; // Time for each iteration
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.4f \t %.4f \n", t, Phi[i]);
        }
        fclose(f);

        /* Save powerspectrum of velocity function to file */
        filename[0] = '\0'; // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/vel_powspec.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            fprintf(f, "%.4f \t %.4f \t %.4f \n", freq[i], Powerspectrum[i], fast_Phi[i]);
        }
        fclose(f);
    }

    /* Free allocated memory */
    free(X);
    X = NULL;
    free(Y);
    Y = NULL;
    free(Z);
    Z = NULL;
    free(Vx);
    Vx = NULL;
    free(Vy);
    Vy = NULL;
    free(Vz);
    Vz = NULL;
    free(Ep);
    Ep = NULL;
    free(Ek);
    Ek = NULL;
    free(Et);
    Et = NULL;
    free(Temperatures);
    Temperatures = NULL;
    free(Pressures);
    Pressures = NULL;
    free(Lat_params);
    Lat_params = NULL;
    free(freq);
    freq = NULL;
    free(Powerspectrum);
    Powerspectrum = NULL;
    free(fast_Phi);
    fast_Phi = NULL;
}

/* Main program */
int main()
{
    /* Initialization of variables and structures */

    double X[N][3];                             // Positions for each particle - x, y, z coordinate.
    int Nc = (int)round(pow(N / 4, 1.0 / 3.0)); // Number of atoms N = 4*Nc*Nc*Nc => Nc = (N/4)^1/3 - 4 atoms per unit cell, with Nc primitive cells in each direction => in total N atoms.
    double m_al = 0.002796;                     // 27 u in our atomic units.
    double dt = 0.005;                          // Recommended according to MD is a few femtoseconds

    /* Task 1 - calculate energy for volumes in the range 64 - 68 Å^3. */
    double v_start = 64;
    double v_end = 68;
    int N_volumes = 9;
    double a_lat = energy_to_volume(N, v_start, v_end, N_volumes, X, Nc);

    /* Task 2 */
    int ndim = 3;
    printf("Equilibrium lattice constant: %.4f Å. \n", a_lat);
    double Teq = 773.15;
    double Peq = 1.0 / 1.602 * 0.000001;
    int equilibrate;
    char label[] = "solid"; // Label for the current production run phase

    /* Code for generating a uniform random number between 0 and 1. srand should only be called once. */
    srand(time(NULL)); // Set the seed for rand
    double rand_disp;  // random displacement

    double x[N][3];         // Current position
    init_fcc(x, Nc, a_lat); // Initialize lattice
    double v[N][3];         // Current velocity
    double a[N][3];         // Current acceleration
    double F[N][3];         // Current force

    /* Generate intial displacements in all directions */
    for (int i = 0; i < N; i -= -1)
    {
        for (int j = 0; j < ndim; j += 1)
        {
            /* Generate random displacement */
            rand_disp = 2 * ((double)rand() / (double)RAND_MAX - 0.5) * 0.065 * a_lat; // Generate random displacement from - a_lat + (+-)0.065*a_lat
            x[i][j] += rand_disp;                                                      // x initialized to be at a_lat
        }
    }

    get_forces_AL(F, x, Nc * a_lat, N); // Calculate initial forces and accelerations
    for (int i = 0; i < N; i -= -1)
    {
        for (int j = 0; j < ndim; j += 1)
        {
            /* Set initial velocities and accelerations */
            v[i][j] = 0.0;
            a[i][j] = F[i][j] / m_al;
        }
    }
    /* Equilibration 1 to melt system */
    //equilibrate = 1;
    //control(x, v, a, F, &a_lat, ndim, Nc, dt, m_al, equilibrate, 1200, Peq, label);
    /* Equilibration 2 to cool down system to 700 K*/
    equilibrate = 1;
    control(x, v, a, F, &a_lat, ndim, Nc, dt, m_al, equilibrate, Teq, Peq, label);

    /* Production */
    equilibrate = 0;
    control(x, v, a, F, &a_lat, ndim, Nc, dt, m_al, equilibrate, Teq, Peq, label);
}
