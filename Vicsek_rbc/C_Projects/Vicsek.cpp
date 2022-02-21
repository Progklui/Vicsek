#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <string>

#include "nrutil.h"
#include "Vicsek.h"

#define PI 3.1415926

Vicsek::Vicsek(double L_init, int N_init, double D_rot_init, double v_init, double R_init, double dt_init, int Nsim_init, int Nsave_init, char * dir_name_init, int init) {
    // constructur: make a system of N spheres in the cubic box
    // x in [-L/2,L/2]
    // y in [-L/2,L/2]
    L        = L_init;

    dt       = dt_init;
    Nsim     = Nsim_init;
    Nsave    = Nsave_init;

    dir_name = dir_name_init;

    // with initial positions on a cubic lattice and L such that the
    // correct packing fraction phi is obtained
    n       = N_init; // int(sqrt((double) N));
    N       = n*n;

    D_rot   = D_rot_init;
    v       = v_init;
    R       = R_init;

    // arrays of length N for positions and directions
    x = dvector(0,N-1);
    y = dvector(0,N-1);

    x_init = dvector(0,N-1);
    y_init = dvector(0,N-1);

    theta = dvector(0,N-1);

    if(init == 1) {
    	printf("\nUniform initialisation.\n\n");
    	init_uniform();
    }

    else {
    	printf("\nRandom initialisation.\n\n");
    	init_random();
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

void Vicsek::rbc(double &x, double &y) {
    // apply periodic boundary conditions
    if (x >= L/2.0) x = -x+L;
    if (x < -L/2.0) x = -x-L;
    if (y >= L/2.0) y = -y+L;
    if (y < -L/2.0) y = -y-L;

    //if (x >= L/2.0) x -= L;
    //if (x < -L/2.0) x += L;
    //if (y >= L/2.0) y -= L;
    //if (y < -L/2.0) y += L;
}

void Vicsek::rbc_theta(double x, double y, double &theta) {
    // apply periodic boundary conditions
    if (x >= L/2.0) theta = PI - theta;
    if (x < -L/2.0) theta = PI - theta;
    if (y >= L/2.0) theta = 2*PI - theta;
    if (y < -L/2.0) theta = 2*PI - theta;
}

double Vicsek::calculate_mean_angle(int i) {
    double avg_sin = 0.; // theta[i]; // 0.;
    double avg_cos = 0.;
    int number_of_particles = 0;

    for (int j=0; j<N; j++) {
        double dx = x[i] - x[j];
        double dy = y[i] - y[j];

        // rbc(dx, dy);

        double r_ij_2 = dx*dx + dy*dy;
        if (r_ij_2 <= R*R) { // && phi < 0.5*phi_vision) {
            // rbc_theta(x[j], y[j], theta[j]);
            avg_sin += sin(theta[j]);
            avg_cos += cos(theta[j]);
            number_of_particles += 1;
        }
    }
    avg_sin = avg_sin/(double) number_of_particles;
    avg_cos = avg_cos/(double) number_of_particles;

    double avg_angle = atan2(avg_sin,avg_cos); // avg_angle/(double) number_of_particles;
    return avg_angle;
}

double Vicsek::check_vision(int i, int j) {
    double dx = x[i] - x[j];
    double dy = y[i] - y[j];

    // rbc(dx, dy);

    double r_ij = sqrt(dx*dx + dy*dy);
    double scalar_prod = dx*cos(theta[i]) + dy*sin(theta[i]);

    double phi = acos(scalar_prod/r_ij);

    return phi;
}

void Vicsek::md_step_vicsek(double dt){
    // update the positions of all particles
    double diff_term = sqrt(2.*D_rot*dt);
    double * avg_angle = dvector(0, N-1);

    for (int i=0; i<N; i++) {
        // calculate average orientation of nearby particles
        avg_angle[i] = calculate_mean_angle(i);
        // avg_angle[i] = check_angle(avg_angle[i]);
    }
    for (int i=0; i<N; i++) {
        // generate noise parameter
        double u1 = rand()/((double) RAND_MAX);
        double u2 = rand()/((double) RAND_MAX);
        double z1 = sqrt(-2.*log(u1))*cos(2.*PI*u2);
        // double z2 = sqrt(-2.*log(u1))*sin(2.*PI*u2);

        //update positions
        x[i] = x[i] + v*cos(theta[i])*dt;
        y[i] = y[i] + v*sin(theta[i])*dt;

        theta[i] = avg_angle[i] + diff_term*z1;
        // theta[i] = fmod(2.*PI + fmod(theta[i], 2.*PI), 2.*PI);
        theta[i] = check_angle(theta[i]);  // check_angle(theta[i], -PI, PI);

        rbc_theta(x[i], y[i], theta[i]);
        theta[i] = check_angle(theta[i]); // atan2(sin(theta[i]), cos(theta[i])); // theta[i] = check_angle(theta[i], -PI, PI);
        // theta[i] = fmod(2.*PI + fmod(theta[i], 2.*PI), 2.*PI);
        rbc(x[i], y[i]);
    }
}

double Vicsek::check_angle(double x) {
    return atan2(sin(x), cos(x)); // x-ceil(x/(2.*PI))*2.*PI+2.*PI; // atan2(sin(x), cos(x));
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

void Vicsek::run_simulation(bool b) {
    // start of the simulation
    char * file_name = "/2_simulation.out";

    printf(" Start of simulation! \n");

    FILE * out;
    std::string buf(dir_name);
    buf.append(file_name);
    out = fopen(buf.c_str(), "w");

    if(b == false) {
    	for (int i=0; i <= Nsim; i++) {
        if (i%Nsave == 0) {
            double va = calculate_va();
            printf(" t = %f, va = %f \n", i*dt, va);
            fprintf(out, "%f %f\n", i*dt, va);
        }
        md_step_vicsek(dt);
    	}
    }
    else {
	    for (int i=0; i <= Nsim; i++) {
		      if (i%Nsave == 0) {
		          double va = calculate_va();
		          printf(" t = %f, va = %f \n", i*dt, va);
		          fprintf(out, "%f %f\n", i*dt, va);
		          store_configuration(i*dt);
		      }
        md_step_vicsek(dt);
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

void Vicsek::init_uniform() {
    // set initial positions and calculate the center of mass
    double cmx = 0.;
    double cmy = 0.;

    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
	          int index = i+j*n;

            x[index] = i;
	          cmx += x[index];

            y[index] = j;
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
    for (int index = 0; index<N; index++) {
        x[index] = (x[index] * L) / ((double) n - 1.);
        y[index] = (y[index] * L) / ((double) n - 1.);

        double u1 = 2*PI*(rand()/((double) RAND_MAX)); // rand()/((double) 2*PI);
        theta[index] = u1;
    }
}

void Vicsek::init_random() {
    // set initial positions and calculate the center of mass
    double cmx = 0.;
    double cmy = 0.;

    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
	          int index = i+j*n;

            x[index] = L * (rand()/((double) RAND_MAX));
	          cmx += x[index];

            y[index] = L * (rand()/((double) RAND_MAX));
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
    for (int index = 0; index<N; index++){
        // x[index] = (x[index] * L) / n;
        // y[index] = (y[index] * L) / n;

        double u1 = 2*PI*(rand()/((double) RAND_MAX)); // rand()/((double) 2*PI);
        theta[index] = u1;
    }
}
