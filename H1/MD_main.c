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


void energy_to_volume(double v_start, double v_end, int N_points, double X[][3], int Nc){
    // TASK 1 - Calculates the energy of N_points volumes between v_start and v_end (INCLUDING end points) 
    // and saves them to a file
    double a_start = pow(v_start, 1.0/3.0);
    double a_end = pow(v_end, 1.0/3.0);
    double energies[N_points];
    double volumes[N_points];
    double delta_a = (a_end-a_start)/((double)N_points-1);
    for (int i=0; i<N_points; i++){
        double a = a_start + delta_a*i;
        init_fcc(X, Nc, a);
        // The supercell is the cell S that describes the same crystal as the unit cell U, but is of larger volume than U.
        // Here, it is the size of the crystal (L = Nc*a).
        double energy = get_energy_AL(X, Nc*a, N);
        // This is wrong - it should be almost much smaller. Are there more unit cells than what I wrote?
        energies[i] = energy/(N/4);  // Energy per unit cell - there are 4 atoms per unit cell in an FCC, and N atoms in total.
        volumes[i] = pow(a, 3.0);
    }  
    // Write energies to file
    FILE *f;
    f = fopen("datafiles/energy_to_volume.dat", "w");
    for(int i=0; i<N_points; i++){
        fprintf(f, "%.4f \t %.4f \n", volumes[i], energies[i]);
    }
    fclose(f);
}


/* Main program */
int main()
{
    // Initialization of variables and structures
    double a0 = 4;  // Lattice parameter in Å
    double X[N][3];  // Positions for each particle - x, y, z coordinate.
    int Nc = (int)round(pow(N/4, 1.0/3.0));  // Number of atoms N = 4*Nc*Nc*Nc => Nc = (N/4)^1/3 - 4 atoms per unit cell, with Nc primitive cells in each direction => in total N atoms. 
    // Initialize lattic
    init_fcc(X, Nc, a0);

    // Code for generating a uniform random number between 0 and 1. srand should only be called once.
    // srand(time(NULL));  // Set the seed for rand
    // double random_value;
    // random_value = (double) rand() / (double) RAND_MAX;
    
    
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
    energy_to_volume(v_start, v_end, N_volumes, X, Nc);
    
    
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
