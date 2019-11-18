/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#define N 256  // Number of atoms
#define timesteps 1000  // Number of timesteps for velocity Verlet


double energy_to_volume(double v_start, double v_end, int N_points, double X[][3], int Nc){
    // TASK 1 - Calculates the energy of N_points volumes between v_start and v_end (INCLUDING end points) 
    // and saves them to a file
    double a_start = pow(v_start, 1.0/3.0);
    double a_end = pow(v_end, 1.0/3.0);
    double energies[N_points];
    double volumes[N_points];
    double delta_a = (a_end-a_start)/((double)N_points-1);
    double lowest_energy = 0;  // Lowest energy 
    double min_a = 0;  // The lattice constant for the minimum energy
    for (int i=0; i<N_points; i++){
        double a = a_start + delta_a*i;
        init_fcc(X, Nc, a);
        // The supercell is the cell S that describes the same crystal as the unit cell U, but is of larger volume than U.
        // Here, it is the size of the crystal (L = Nc*a).
        double energy = get_energy_AL(X, Nc*a, N);
        // This is wrong - it should be almost much smaller. Are there more unit cells than what I wrote?
        energies[i] = energy/(N/4);  // Energy per unit cell - there are 4 atoms per unit cell in an FCC, and N atoms in total.
        volumes[i] = pow(a, 3.0);
        if(energies[i] < lowest_energy){
            lowest_energy = energies[i];
            min_a = a;
        } 
    }  
    // Write energies to file
    FILE *f;
    f = fopen("datafiles/energy_to_volume.dat", "w");
    for(int i=0; i<N_points; i++){
        fprintf(f, "%.4f \t %.4f \n", volumes[i], energies[i]);
    }
    fclose(f);
    return(min_a); 
}


void velocity_verlet(double a0, int ndim, int Nc, double dt, double m_al){
    // Lattice constant a, number of dimensions, number of unit cells in all directions and mass of Al. 
    // Modified for task 5 to calculate the mean squared displacement as a function of time

    // Code for generating a uniform random number between 0 and 1. srand should only be called once.
    srand(time(NULL));  // Set the seed for rand
    double rand_disp;  // random displacement
    

    // Initialization of variables
    int i,j,k; // loop variables 
    double x[N][3];  // Current position
    init_fcc(x, Nc, a0);  // Initialize lattice
    double v[N][3];  // Current velocity
    double a[N][3];  // Current acceleration
    double F[N][3];  // Current force

    double Ep[timesteps+1];  // Potential energy
    double Ek[timesteps+1];  // Kinetic energy
    double Et[timesteps+1];  // Total energy

    double MSD[timesteps+1];  // Task 5 - Mean squared displacement as measured from the start point
    double phi[timesteps+1];  // Task 6 - Time autocorrelation function - delta t is ten iterations

    // Task 5 - mean squared displacement
    double x0[N][3];
    for(i=0; i<N; i++){
        for(j=0; j<ndim; j++){
            x0[i][j] = x[i][j];  // x0 is the reference positions for all atoms
        }
    } 
    
    // Generate intial displacements in all directions
    for(i=0; i<N; i-=-1){
        for(j=0; j<ndim; j+=1){
            // Generate random displacement
            rand_disp = 2*((double) rand() / (double) RAND_MAX  - 0.5)*0.065*a0;  // Generate random displacement from - a0 + (+-)0.065*a0
            x[i][j] +=  rand_disp;  // x initialized to be at a0
        }
        
    }

    get_forces_AL(F, x, Nc*a0, N);  // Calculate initial forces and accelerations
    for(i=0; i<N+1; i-=-1){
        for(j=0; j<ndim; j+=1){
            // Set initial velocities and accelerations
            v[i][j] = 0.0; 
            a[i][j] = F[i][j]/m_al;
        }
    }

    // Check that initial velocities are indeed 0
    for(i=0; i<N; i-=-1){
        for(j=0; j<ndim; j+=1){
            if(v[i][j] != 0.0){
                printf("Initialization error: %.2f \n", v[i][j]);
            }
        }
    }

    Ep[0] = get_energy_AL(x, Nc*a0, N);  // Supercell length is Nc*a0
    // Calculate the kinetic energy - T = 0.5*m*|v|^2
    double ek = 0;
    double abs_v_sq;
    for(j=0; j<N; j+=1){
        abs_v_sq = 0;
        for(k=0; k<ndim; k+=1){
            abs_v_sq += v[j][k] * v[j][k];  // ek = 0.5*mv^2 
        }
        // printf("v[0]: %.2f \t v[1]: %.2f \t  v[2]: %.2f \n", v[j][0], v[j][1], v[j][2]);
        ek += 0.5*m_al*abs_v_sq;
    }
    Ek[0] = ek; 
    Et[0] = Ep[0] + Ek[0];  // Total energy is sum of kinetic and potential
    

    // Velocity verlet 
    for(i=1; i<timesteps+1; i+=1){
        // Calculate approximate new velocity - at deltat/2
        for(j=0; j<N; j+=1){
            for(k=0; k<ndim; k+=1){
                v[j][k] += dt/2 * a[j][k];
            }
        }

        // Calculate new positions
        for(j=0; j<N; j+=1){
            for(k=0; k<ndim; k+=1){
                x[j][k] += dt * v[j][k];
            }
        }

        // Calculate new forces and accelerations
        get_forces_AL(F, x, Nc*a0, N);  // Supercell length is Nc*a0
        for(j=0; j<N; j+=1){
            for(k=0; k<ndim; k+=1){
                a[j][k] = F[j][k]/m_al;
                // printf("a[0]: %.2f \t a[1]: %.2f \t  a[2]: %.2f \n", a[j][0], a[j][1], a[j][2]);
                // a[j][k] = 0.0;
            }
        }

        // Calculate final velocities
        for(j=0; j<N; j+=1){
            for(k=0; k<ndim; k+=1){
                v[j][k] += dt/2 * a[j][k];
            }
        }

        // Calculate and save energies for this iteration
        Ep[i] = get_energy_AL(x, Nc*a0, N);  // Supercell length is Nc*a0
        // Calculate the kinetic energy - T = 0.5*m*|v|^2
        double ek = 0;
        double abs_v_sq;
        for(j=0; j<N; j+=1){
            abs_v_sq = 0;
            for(k=0; k<ndim; k+=1){
                abs_v_sq += v[j][k] * v[j][k];  // ek = 0.5*mv^2 
            }
            // printf("v[0]: %.2f \t v[1]: %.2f \t  v[2]: %.2f \n", v[j][0], v[j][1], v[j][2]);
            ek += 0.5*m_al*abs_v_sq;
        }
        Ek[i] = ek; 
        Et[i] = Ep[i] + Ek[i];  // Total energy is sum of kinetic and potential

        // Task 5 - Calculate MSD for this iteration
        double sum = 0;
        double delta_x;
        for(j=0; j<N; j+=1){
            delta_x = 0;
            for(k=0; k<ndim; k+=1){
                delta_x += pow(x[j][k]-x0[j][k], 2.0);
            }
            // printf("dx[0]: %.2f \t dx[1]: %.2f \t  dx[2]: %.2f \n", x[j][0]-x0[j][0], x[j][1]-x0[j][1], x[j][2]-x0[j][2]);
            sum += delta_x;
        }
        MSD[i] = sum/N;
    }

    // Save energies to file
    FILE *f;
    f = fopen("datafiles/vv_energies.dat", "w");
    double t;  // Time for each iteration
    for(i=0; i<timesteps+1; i++){
        t = i*dt;
        fprintf(f, "%.4f \t %.4f \t %.4f \t %.4f \n", t, Ep[i], Ek[i], Et[i]);
    }
    fclose(f);
    // Saved MSD to file
    f = fopen("datafiles/MSD.dat", "w");
    t=0;  // Time for each iteration
    for(i=0; i<timesteps+1; i++){
        t = i*dt;
        fprintf(f, "%.4f \t %.4f \n", t, MSD[i]);
    }
    fclose(f);

}

/* Main program */
int main()
{
    // Initialization of variables and structures
    
    double X[N][3];  // Positions for each particle - x, y, z coordinate.
    int Nc = (int)round(pow(N/4, 1.0/3.0));  // Number of atoms N = 4*Nc*Nc*Nc => Nc = (N/4)^1/3 - 4 atoms per unit cell, with Nc primitive cells in each direction => in total N atoms. 
    double m_al = 0.002796;  // 27 u in our atomic units.
    double dt = 0.005;  // Recommended according to MD is a few femtoseconds
    // Initialize lattic
    // double a0 = 4;  // Lattice parameter in Å
    // init_fcc(X, Nc, a0);

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

    // Task 1 - calculate energy for volumes in the range 64 - 68 Å^3. 
    double v_start = 64;
    double v_end = 68;
    int N_volumes = 9;
    double min_a;
    min_a = energy_to_volume(v_start, v_end, N_volumes, X, Nc);
    
    // Task 2
    int ndim = 3;
    printf("Equilibrium lattice constant: %.4f Å. \n", min_a);
    velocity_verlet(min_a, ndim, Nc, dt, m_al);
    
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
    
    
    
}
