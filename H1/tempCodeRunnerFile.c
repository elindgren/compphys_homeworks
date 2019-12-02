    // if (!equilibrate)
    // {   
    //     omp_set_num_threads(8);
    //     #pragma omp parallel sections
    //     {
    //         #pragma omp section
    //         {
    //             printf("Calculating MSD function \n");
    //             mean_squared_displacement(timesteps+1, N, MSD, X, Y, Z);  // Task 5 
    //             printf("Finished MSD function \n");
    //         }
    //         #pragma omp section
    //         {
    //             printf("Calculating velocity correlation function \n");
    //             velocity_correlation(timesteps+1, N, Phi, Vx, Vy, Vz); // Task 6
    //             printf("Finished velocity correlation function \n");
    //         }
    //         #pragma omp section
    //         {
    //             printf("Calculating power spectrum function \n");
    //             fast_velocity_correlation(timesteps+1, N, fast_Phi, Powerspectrum, freq, Vx, Vy, Vz, dt);  // Task 7 
    //             printf("Finished power spectrum function \n");
    //         }
    //     }
    // }