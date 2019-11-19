#include <math.h>
#include "alpotential.h"


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

    return m / (3 * N * k ) * sum;
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
    sum *= m/3.0;

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
