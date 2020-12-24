/*
 * main.cpp
 *
 *  Created on: Dec 20, 2020
 *      Author: d-w-h
 *
 *      Numerical solution of the Schrodinger equation:
 *
 *      ih*d(psi)/dt = -h*h/2m*d^2(psi)/dx^2 + V*psi
 */

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "complex.hpp"
#include "solver.hpp"
#include "user_types.hpp"

double V(double x, double t) {
	/* potential */
	double k = 0.1;
	return 0.5*k*x*x;
}

int main(int argc, char* argv[]) {
    d_data domain_data;
    t_data time_data;
    p_params physical_params;
    bound_and_psi_init boundaries_and_psi_init;
    s_data solver_data;
    clock_t t_start = clock();

    /* Parameters */
    domain_data.Nx = 39;                           //Number of nodes along x axis, should be an odd number
    domain_data.Nt = 600;                          //Number of timesteps
    domain_data.L = 15.0;                          //Length of domain

    time_data.to = 0.0;                            //Initial time
    time_data.tf = 30.0;                           //Final time

    physical_params.h = 1.0;                       //Constant
    physical_params.m = 1.0;                       //Mass of particle

    boundaries_and_psi_init.psi_west.a = 0.0;      //Real part of west boundary
    boundaries_and_psi_init.psi_west.b = 0.0;      //Imaginary part of west boundary
    boundaries_and_psi_init.psi_east.a = 0.0;      //Real part of east boundary
    boundaries_and_psi_init.psi_east.b = 0.0;      //Imaginary part of east boundary
    boundaries_and_psi_init.psi_init.a = 0.1;      //Real value of the wave function at t = 0
    boundaries_and_psi_init.psi_init.b = 0.5;      //Imaginary value of the wave function at t = 0

    /* Allocate memory for solver results */
    solver_data.psi = new Complex[domain_data.Nx];
    solver_data.x_c = new double[domain_data.Nx];

    /* Execute solver */
    solver(domain_data,
           time_data,
           physical_params,
           boundaries_and_psi_init,
           &solver_data);

    /* Print results */
    for(int j = 0; j < domain_data.Nx; ++j) {
        printf("x_c[%i]: %f, real(psi[%i]): %f, Im(psi[%i]): %f\n",
                j, solver_data.x_c[j], j, solver_data.psi[j].a, j, solver_data.psi[j].b);
    }

    /* Deallocate memory for solver results */
    delete [] solver_data.psi;
    delete [] solver_data.x_c;

    /* Print time taken for execution */
    clock_t t_end = clock();
    double sim_time = double (t_end - t_start)/CLOCKS_PER_SEC;
    printf("time taken: %f\n", sim_time);

    return 0;
}
