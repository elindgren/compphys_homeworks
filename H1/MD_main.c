/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
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
#define N 256          // Number of atoms
#define timesteps 5000 // Number of timesteps for velocity Verlet

/*
 * Encapsulates the velocity Verlet algorithm
 * Lattice constant a, number of dimensions, number of unit cells in all directions and mass of Al
 * Modified for task 5 to calculate the mean squared displacement as a function of time
 */
void velocity_verlet(double x[][3], double v[][3], double a[][3], double F[][3], double *a_lat, int ndim, int Nc, double dt, double m_al, int equilibrate, double Teq, double Peq, char label[])
{
    int i, j, k; // loop variables

    double Ep[timesteps + 1]; // Potential energy
    double Ek[timesteps + 1]; // Kinetic energy
    double Et[timesteps + 1]; // Total energy

    /* Task 3 */
    double T, P;                 // Pressure and temperature
    double alpha_t, alpha_p;     // Scaling parameters for equilibration
    double timedecay = 250 * dt; // Timedecay
    double Temperatures[timesteps + 1];
    double Pressures[timesteps + 1];
    double Lat_params[timesteps + 1];
    double positions[timesteps + 1][ndim];

    /* Task 5 */
    double MSD[timesteps + 1]; // Mean squared displacement as measured from the start point
    
    double *x_pos[timesteps + 1];   // Matrix containing all particle x positions
    for (i = 0; i < timesteps + 1; i++){
        x_pos[i] = (double *)malloc(N * sizeof(double));
        for (j = 0; j < N; j++){
            x_pos[i][j] = 0.0;
        }
    }

    double x0[N][3];           // Save reference positions in x0
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            x0[i][j] = x[i][j]; // x0 is the reference positions for all atoms
        }
    }

    /* Task 6 */
    double phi[timesteps + 1];                      // Time velocity correlation function - not defined for last corr_offset values.
    double *x_vel[timesteps + 1];   // Matrix containing all particle x velocities
    for (i = 0; i < timesteps + 1; i++){
        x_vel[i] = (double *)malloc(N * sizeof(double));
        for (j = 0; j < N; j++){
            x_vel[i][j] = 0.0;
        }
    }

    /* Task 7 */
    int B = 2*(timesteps+1);  // N in the notation in the notebook
    double h[B];  // Vector to hold raw data for a single particle
    double H[B];  // FFT of one particle
    double C[timesteps+1];  // Correlation function from inverse FFT
    double freq[B];
    double *x_vel_pad[B];  // Same as x_vel, but padded with N+1 zeros at the end
    double *H_all[B];  // Fourier transform of velocity for all particles
    for (i = 0; i < B; i++){
        x_vel_pad[i] = (double *)malloc(N * sizeof(double));
        H_all[i] = (double *)malloc(N * sizeof(double));
        for (j = 0; j < N; j++){
            x_vel_pad[i][j] = 0.0;
            H_all[i][j] = 0.0;
        }
    }

    /* Calculate energies for initial conditions */
    Ep[0] = get_energy_AL(x, Nc * *a_lat, N); // Supercell length is Nc*a_lat

    /* Calculate the kinetic energy - T = 0.5*m*|v|^2 */
    double ek = 0;
    double abs_v_sq;
    for (j = 0; j < N; j += 1)
    {
        abs_v_sq = 0;
        for (k = 0; k < ndim; k += 1)
        {
            abs_v_sq += v[j][k] * v[j][k]; // ek = 0.5*mv^2
        }
        // printf("v[0]: %.2f \t v[1]: %.2f \t  v[2]: %.2f \n", v[j][0], v[j][1], v[j][2]);
        ek += 0.5 * m_al * abs_v_sq;
    }
    Ek[0] = ek;
    Et[0] = Ep[0] + Ek[0]; // Total energy is sum of kinetic and potential

    /* Save initial postion for particle p */
    int p = 0;
    for (j = 0; j < ndim; j++)
    {
        positions[0][j] = x[0][j];
    }
    /* Save initial positions and velocities for particles */
    for (j = 0; j < N; j += 1)
    {   
        x_pos[0][j] = x[j][0];
        x_vel[0][j] = v[j][0]; // We only take x direction
    }

    /* The velocity Verlet algorithm */
    for (i = 1; i < timesteps + 1; i += 1)
    {
        /* v(t+dt/2) */
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < ndim; k += 1)
            {
                v[j][k] += dt / 2 * a[j][k];
            }
        }

        /* q(t+dt) */
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < ndim; k += 1)
            {
                x[j][k] += dt * v[j][k];
            }
        }

        /* a(t+dt) */
        get_forces_AL(F, x, Nc * *a_lat, N); // Supercell length is Nc*a_lat
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < ndim; k += 1)
            {
                a[j][k] = F[j][k] / m_al;
                // printf("a[0]: %.2f \t a[1]: %.2f \t  a[2]: %.2f \n", a[j][0], a[j][1], a[j][2]);
                // a[j][k] = 0.0;
            }
        }

        /* v(t+dt) */
        for (j = 0; j < N; j += 1)
        {
            for (k = 0; k < ndim; k += 1)
            {
                v[j][k] += dt / 2 * a[j][k];
            }
        }

        /* Save the position of particle p */
        for (j = 0; j < ndim; j++)
        {
            positions[i][j] = x[p][j];
        }

        /* Calculate and save energies for this iteration */
        Ep[i] = get_energy_AL(x, Nc * *a_lat, N); // Supercell length is Nc*a_lat
        /* Calculate the kinetic energy - T = 0.5*m*|v|^2 */
        ek = 0;
        for (j = 0; j < N; j += 1)
        {
            abs_v_sq = 0;
            for (k = 0; k < ndim; k += 1)
            {
                abs_v_sq += v[j][k] * v[j][k]; // ek = 0.5*mv^2
            }
            // printf("v[0]: %.2f \t v[1]: %.2f \t  v[2]: %.2f \n", v[j][0], v[j][1], v[j][2]);
            ek += 0.5 * m_al * abs_v_sq;
        }
        Ek[i] = ek;
        Et[i] = Ep[i] + Ek[i]; // Total energy is sum of kinetic and potential

        /* Task 5 - Calculate MSD for this iteration */
        // TODO change to same format as absolute velocities - save absolute displacement squared
        double sum = 0;
        double delta_x;
        for (j = 0; j < N; j += 1)
        {
            delta_x = 0;
            for (k = 0; k < ndim; k += 1)
            {
                delta_x += pow(x[j][k] - x0[j][k], 2.0);
            }
            // printf("dx[0]: %.2f \t dx[1]: %.2f \t  dx[2]: %.2f \n", x[j][0]-x0[j][0], x[j][1]-x0[j][1], x[j][2]-x0[j][2]);
            sum += delta_x;
        }
        MSD[i] = sum / N;

        /* Task 6 - Calculate velocity correlation function */
        /* Save current positions and velocities for autocorrelation functions */
        for (j = 0; j < N; j += 1)
        {
            x_pos[i][j] = x[j][0];
            x_vel[i][j] = v[j][0];
            x_vel_pad[i][j] = v[j][0];
        }

        // Set scaling parameters
        T = calc_temp(N, m_al, v);
        Temperatures[i] = T;
        P = calc_pres(N, m_al, v, x, Nc * *a_lat);
        Pressures[i] = P;
        if (equilibrate)
        {
            alpha_t = calc_temp_scale(dt, timedecay, Teq, T);
            alpha_p = calc_pres_scale(dt, timedecay, Peq, P);
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
    if (!equilibrate){
        // printf("Calculating MSD function \n");
        /* I tried to do the same with MSD as with velocity correlation */
        // calc_corr_function(timesteps + 1, N, MSD, x_pos);  // Task 5
        printf("Calculating velocity correlation function \n");
        calc_corr_function(timesteps + 1, N, phi, x_vel);  // Task 6
    }

    /* Task 7 - Fourier Transform velocities */
    // Using FFT-module from E1 
    if (!equilibrate){
        fft_freq(freq, dt, B);  // Corresponding frequencies to powerspectra
        for(i=0; i<N; i++){
            /* Calculate power spectrum for one particle */
            for(j=0; j<B; j++){
                h[j] = x_vel_pad[j][i];  // Extract velocities for a single particle
            }
            powerspectrum(h, H, B);  // Powerspectrum
            // Write powerspectrum to H_all
            for(j=0; j<B; j++){
                H_all[j][i] = H[j];  // Write spectrum for particle i
            }
            /* Perform inverse FFT */
            for(int l=0; l<timesteps+1; l++){
                double Cl = 0;
                for(int m=0; m<N; m++){
                    // Cl
                }
            }
        }

    }


    FILE *f;
    double t;
    char filename[50] = "";
    char dir[100] = "datafiles/";
    strcat(dir, label);
    strcat(filename, dir);
    if (equilibrate)
    {   
        /* Equilibration run */
        printf("Saving equilibration \n");

        filename[0] = '\0';  // Empty filename string
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
        filename[0] = '\0';  // Empty filename string
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
        filename[0] = '\0';  // Empty filename string
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
        filename[0] = '\0';  // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/MSD.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", t, positions[i][0], positions[i][1], positions[i][2], MSD[i]);
        }
        fclose(f);

        /* Save velocity correlation function to file */
        filename[0] = '\0';  // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/vel_corr.dat");
        f = fopen(filename, "w");
        t = 0; // Time for each iteration
        for (i = 0; i < timesteps + 1; i++)
        {
            t = i * dt;
            fprintf(f, "%.4f \t %.4f \n", t, phi[i]);
        }
        fclose(f);

        /* Save powerspectrum of velocity function to file */
        filename[0] = '\0';  // Empty filename string
        strcat(filename, dir);
        strcat(filename, "/vel_powspec.dat");
        f = fopen(filename, "w");
        for (i = 0; i < timesteps + 1; i++)
        {
            fprintf(f, "%.4f ", freq[i]);
            for(j=0; j<N; j++){
                fprintf(f, "\t %.4f ", H_all[i][j]);
            }
            fprintf(f, " \n ");
        }
        fclose(f);
    }
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
    char label[] = "solid";  // Label for the current production run phase

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
    // equilibrate = 1;
    // velocity_verlet(x, v, a, F, &a_lat, ndim, Nc, dt, m_al, equilibrate, 1200, Peq, label);
    /* Equilibration 2 to cool down system to 700 K*/
    equilibrate = 1;  
    velocity_verlet(x, v, a, F, &a_lat, ndim, Nc, dt, m_al, equilibrate, Teq, Peq, label);

    /* Production */
    equilibrate = 0;
    velocity_verlet(x, v, a, F, &a_lat, ndim, Nc, dt, m_al, equilibrate, Teq, Peq, label);
}

/*
 * OLD FRICKING COMMENTS!
 */

/*
     Descriptions of the different functions in the files initfcc.c and
     alpotential.c are listed below.
    */

/* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */

/* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
/*
     double energy;
     energy = get_energy_AL(pos, L, N);
    */

/* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
/*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */

/*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
/*
     get_forces_AL(f,pos, L, N);
    */
