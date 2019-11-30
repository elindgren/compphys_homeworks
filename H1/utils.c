#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "alpotential.h"
#include "initfcc.h"
#include "fft_func.h"

/* Calculates the potential energy as a function of the cell volume.
 * Saves energy to file "energy_to_volume.dat"
 */
double energy_to_volume(int N, double v_start, double v_end, int N_points, double X[][3], int Nc)
{
    double a_start = pow(v_start, 1.0 / 3.0);
    double a_end = pow(v_end, 1.0 / 3.0);
    double energies[N_points];
    double volumes[N_points];
    double delta_a = (a_end - a_start) / ((double)N_points - 1);
    double lowest_energy = 0; // Lowest energy
    double min_a = 0;         // The lattice constant for the minimum energy
    double a;
    double energy;
    FILE *f;

    for (int i = 0; i < N_points; i++)
    {
        a = a_start + delta_a * i;
        init_fcc(X, Nc, a);
        /* The supercell is the cell S that describes the same crystal 
         * as the unit cell U, but is of larger volume than U.
         * Here, it is the size of the crystal (L = Nc*a).
         */
        energy = get_energy_AL(X, Nc * a, N);
        energies[i] = energy / (N / 4); /* Energy per unit cell - there are 4 atoms 
                                         * per unit cell in an FCC, and N atoms in total.
                                         */
        volumes[i] = pow(a, 3.0);
        if (energies[i] < lowest_energy)
        {
            lowest_energy = energies[i];
            min_a = a;
        }
    }

    /* Write energies to file */
    f = fopen("datafiles/energy_to_volume.dat", "w");
    for (int i = 0; i < N_points; i++)
    {
        fprintf(f, "%.4f \t %.4f \n", volumes[i], energies[i]);
    }
    fclose(f);
    return (min_a);
}

/* 
 * Calculate a general correlation function.
 * Phi is the output correlation function at various time t.
 * A is a matrix containing the observable in question for all particles at various times t.
 * M is the total number of timesteps.
 * corr_offset is the offset at which the correlation function is to be calculated.
 */
void calc_corr_function(int M, int N, double phi[], double A[][N])
{
    // Initialize variables
    int l, m, n; // Iteration variables

    // double corr_func; // The correlation at time t for one particle
    double sum; // Sum of correlations at time t for all particles
    double Cn;  // correlation function value at time t for particle n
    for (l = 0; l < M; l += 1)
    {   
        // printf("\t l = %d \n", l);
        // l is timelag in number of timesteps
        sum = 0;
        // For each timelag l, calculate the correlation function for this timestep for all particles n
        for (n = 0; n < N; n++){
            Cn = 0;
            for (m = 0; m < M-l; m += 1)
            {
                // printf("\t m = %d \n", m);
                Cn += A[m+l][n] * A[m][n];  // Particle 0 atm
            }
            // Average over times M-l
            sum += Cn / (M - l);
        }
        /* Average over all particles */
        phi[l] = sum / N;
    }
}

void mean_squared_displacement(int M, int N, double delta[], double X[][N], double Y[][N], double Z[][N])
{
    int l, m, n; // Iteration variables

    // double corr_func; // The correlation at time t for one particle
    double sum; // Sum of correlations at time t for all particles
    double delta_n;  // correlation function value at time t for particle n
    for (l = 0; l < M; l += 1)
    {   
        /* l is timelag in number of timesteps */
        if(l%1000==0){
            printf("\tMSD: %.2f \n", (double)l/M);
        }
        sum = 0;
        /* For each timelag l, calculate the correlation function for this timestep for all particles n */
        for (n = 0; n < N; n++){
            delta_n = 0;
            for (m = 0; m < M-l; m += 1)
            {
                delta_n += pow(X[m+l][n] - X[m][n], 2.0) + pow(Y[m+l][n] - Y[m][n], 2.0) + pow(Z[m+l][n] - Z[m][n], 2.0);  // Particle 0 atm
            }
            /* Average over times M-l */
            sum += delta_n / (M - l);
        }
        /* Average over all particles */
        delta[l] = sum / N;
    }
}

void velocity_correlation(int M, int N, double phi[], double Vx[][N], double Vy[][N], double Vz[][N])
{
    int l, m, n; // Iteration variables

    // double corr_func; // The correlation at time t for one particle
    double sum; // Sum of correlations at time t for all particles
    double phi_n;  // correlation function value at time t for particle n
    for (l = 0; l < M; l += 1)
    {   
        /* l is timelag in number of timesteps */
        if(l%1000==0){
            printf("\tVelocity corr: %.2f \n", (double)l/M);
        }
        sum = 0;
        /* For each timelag l, calculate the correlation function for this timestep for all particles n */
        for (n = 0; n < N; n++){
            phi_n = 0;
            for (m = 0; m < M-l; m += 1)
            {
                phi_n += Vx[m+l][n]*Vx[m][n] + Vy[m+l][n]*Vy[m][n] + Vz[m+l][n]*Vz[m][n];  // Particle 0 atm
            }
            /* Average over times M-l */
            sum += phi_n / (M - l);
        }
        /* Average over all particles */
        phi[l] = sum / N;
    }
}

void fast_velocity_correlation(int timesteps, int N, double fast_phi[], double P[], double freq[], double Vx[][N], double Vy[][N], double Vz[][N], double dt){
    int B = 2*timesteps;
    // Declare matrices for FFT
    double *vx = malloc(B * sizeof(double));  // Save one velocity vector
    double *vy = malloc(B * sizeof(double));  // Save one velocity vector
    double *vz = malloc(B * sizeof(double));  // Save one velocity vector

    double *Px = malloc(B * sizeof(double));  // Corresponding FFT vector
    double *Py = malloc(B * sizeof(double));  // Corresponding FFT vector
    double *Pz = malloc(B * sizeof(double));  // Corresponding FFT vector


    fft_freq(freq, dt, B);  // Corresponding frequencies to powerspectra
    for(int i = 0; i<N; i++){
        /* Save to vector v */
        for(int j = 0; j<B; j++){
            if(j<timesteps){
                vx[j] = Vx[j][i];
                vy[j] = Vy[j][i];
                vz[j] = Vz[j][i];
            }else{
                /* Pad with zeros */
                vx[j] = 0.0;
                vy[j] = 0.0;
                vz[j] = 0.0;
            }
        }
        /* Calculate FFT */
        powerspectrum(vx, Px, B);
        powerspectrum(vy, Pz, B);
        powerspectrum(vz, Pz, B);
        /* Save to H matrix */
        for(int j = 0; j<B; j++){
            P[j] += Px[j] + Py[j] + Pz[j];
        }
    }

    /* Calculate average */
    for(int j = 0; j<B; j++){
        P[j] /= (double)N;
    }
    /* Calculate powerspectrum of Powerspectrum */
    fft(P, fast_phi, B);
    /* Square root since we want fourier transform of power spectrum, not powerspectrum. */
    for(int j = 0; j<B; j++){
        fast_phi[j] = fast_phi[j];
    }

    free(vx); vx=NULL; free(vy); vy=NULL; free(vz); vz=NULL;
    free(Px); Px=NULL; free(Px); Py=NULL; free(Pz); Pz=NULL;
}


/* Calculate the instantaneous temperature based on kinetic energy */
double calc_temp(int N, double m, double v[][3])
{
    double k = 0.000086169; // Boltzmann constant in atomic units
    double sum = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum += v[i][j] * v[i][j];
        }
    }

    return m / (3 * N * k) * sum;
}

/* Calculate the instantaneous pressure bassed on the kinetict energy and virial */
double calc_pres(int N, double m, double v[][3], double x[][3], double cell_length)
{
    double sum;
    double volume;

    /* Sum the contribution from kinetic energy */
    sum = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum += v[i][j] * v[i][j];
        }
    }
    sum *= m / 3.0;

    /* Add the virial contribution */
    volume = pow(cell_length, 3.0); // Volume of supercell
    sum += get_virial_AL(x, cell_length, N);
    sum /= volume;

    return sum;
}

/* Calculate the velocity scaling parameter for the temperature scaling */
double calc_temp_scale(double dt, double time_decay, double temp_eq, double current_temp)
{
    return 1 + 2 * dt / time_decay * (temp_eq - current_temp) / current_temp;
}

/* Calculate the position and box size scaling parameter for the pressure scaling */
double calc_pres_scale(double dt, double time_deacy, double pres_eq, double current_pres)
{
    double kappa = 160.2 / 79.0; // In units Å^3 / eV
    return 1 - kappa * dt / time_deacy * (pres_eq - current_pres);
}
