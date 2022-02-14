#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <string>

#include "nrutil.h"
#include "Vicsek.h"

#define PI 3.1415926

Vicsek::Vicsek(double sigma_init, int N_init, double phi, double D_rot_init, double v_init, double R_init, double dt_init, int Nsim_init, int Nsave_init, char * dir_name_init) {
    // constructur: make a system of N spheres in the cubic box
    // x in [-L/2,L/2]
    // y in [-L/2,L/2
    dt       = dt_init;
    Nsim     = Nsim_init;
    Nsave    = Nsave_init;

    dir_name = dir_name_init;

    // with initial positions on a cubic lattice and L such that the
    // correct packing fraction phi is obtained
    N       = N_init;
    n       = int (pow(N+0.000002,1./3));
    sigma   = sigma_init;
    sigma2  = sigma*sigma;

    D_rot   = D_rot_init;
    v       = v_init;
    R       = R_init;

    // arrays of length N for positions and directions
    x = dvector(0,N-1);
    y = dvector(0,N-1);

    x_init = dvector(0,N-1);
    y_init = dvector(0,N-1);

    theta = dvector(0,N-1);

    // set initial positions and calculate the center of mass
    double cmx = 0.;
    double cmy = 0.;

    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
	          int index = i+j*n;

            x[index] = i*sigma;
	          cmx += x[index];

            y[index] = j*sigma;
	          cmy += y[index];
	      }
    }
    cmx = cmx/N;
    cmy = cmy/N;

    for (int index = 0; index<N; index++) {
        x[index] -= cmx;
        y[index] -= cmy;
    }

    // rescale the system to match the packing fraction
    double Lo = n*sigma;
    L = pow(N*PI*sigma*sigma*sigma/6/phi,1./3);
    for (int index = 0; index<N; index++){
        x[index] = x[index]*L/Lo;
        y[index] = y[index]*L/Lo;

        double u1 = 2*PI*(rand()/((double) RAND_MAX)); // rand()/((double) 2*PI);
        theta[index] = u1;
    }

    save_simulation_params();
}

void Vicsek::pbc(double &x, double &y) {
    // apply periodic boundary conditions
    if (x >= L/2.0) x -= L;
    if (x < -L/2.0) x += L;
    if (y >= L/2.0) y -= L;
    if (y < -L/2.0) y += L;
}

double Vicsek::get_size() {
    // return the size of the system
    return L;
}

double Vicsek::get_packing_fraction() {
    // return the size of the system
    return N*PI*sigma*sigma*sigma/6./(L*L);
}

double Vicsek::calculate_mean_angle(int i) {
    double avg_angle = 0.;
    int number_of_particles = 0;

    for (int j=0; j<N; j++) {
        double r_ij_2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
        if (r_ij_2 < R*R) {
            avg_angle += theta[j];
            number_of_particles += 1;
        }
    }
    avg_angle = avg_angle/number_of_particles;
    return avg_angle;
}

void Vicsek::md_step_vicsek(double dt){
    // update the positions of all particles
    double diff_term = sqrt(2.*D_rot*dt);
    for (int i=0; i<N; i++) {
        // generate noise parameter
        double u1 = rand()/((double) RAND_MAX);
        double u2 = rand()/((double) RAND_MAX);
        double z1 = sqrt(-2.*log(u1))*cos(2.*PI*u2);
        double z2 = sqrt(-2.*log(u1))*sin(2.*PI*u2);

        // calculate average orientation of nearby particles
        double avg_angle = calculate_mean_angle(i);

        //update positions
        x[i] = x[i] + v*cos(theta[i])*dt;
        y[i] = y[i] + v*sin(theta[i])*dt;

        theta[i] = avg_angle + diff_term*z1;

        pbc(x[i], y[i]);
    }
}

double Vicsek::calculate_va() {
    // calculate normalized average velocity
    double vx = 0.;
    double vy = 0.;

    for (int i=0; i<N; i++) {
        vx += cos(theta[i]);
        vy += sin(theta[i]);
    }
    double va = sqrt(vx*vx + vy*vy) / N;
    return va;
}

void Vicsek::md_equilibration(int Neq) {
    // equilibrate the system
    int nsave = Neq/100;
    char * file_name = "/1_equilibration.out";

    printf(" Start of equilibration! \n");

    FILE * out;
    std::string buf(dir_name);
    buf.append(file_name);
    out = fopen(buf.c_str(), "w");

    for (int i=1; i <= Neq; i++) {
        md_step_vicsek(dt);
        if (i%nsave == 0) {
            double va = calculate_va();
            printf(" t = %f, va = %f \n", i*dt, va);
            fprintf(out, "%f %f\n", i*dt, va);
        }
    }
    fclose(out);
    printf(" Equilibration done! \n\n");
}

void Vicsek::store_configuration(double t) {
    std::string file_name = "/configuration_t_" + std::to_string(t);

    FILE * out;
    std::string buf(dir_name);
    buf.append(file_name);
    out = fopen(buf.c_str(), "w");
    for (int i=0; i<N; i++) {
        fprintf(out, "%f %f %f\n", x[i], y[i], theta[i]);
    }
    fclose(out);
}

void Vicsek::run_simulation() {
    // start of the simulation
    char * file_name = "/2_simulation.out";

    printf(" Start of simulation! \n");

    FILE * out;
    std::string buf(dir_name);
    buf.append(file_name);
    out = fopen(buf.c_str(), "w");

    for (int i=1; i <= Nsim; i++) {
        md_step_vicsek(dt);
        if (i%Nsave == 0) {
            double va = calculate_va();
            printf(" t = %f, va = %f \n", i*dt, va);
            fprintf(out, "%f %f\n", i*dt, va);
            store_configuration(i*dt);
        }
    }
    fclose(out);
    printf(" Simulation done! \n\n");
}

void Vicsek::save_simulation_params() {
    char * file_name = "/0_params.out";

    FILE * out;
    std::string buf(dir_name);
    buf.append(file_name);
    out = fopen(buf.c_str(), "w");

    fprintf(out, "N, %i\n", N);
    fprintf(out, "L, %f\n", L);
    fprintf(out, "dt, %f\n", dt);
    fprintf(out, "R, %f\n", R);
    fprintf(out, "v, %f\n", v);
    fprintf(out, "D, %f\n", D_rot);
    fprintf(out, "Nsim, %i\n", Nsim);
    fprintf(out, "Nsav, %i\n", Nsave);

    fclose(out);
}
