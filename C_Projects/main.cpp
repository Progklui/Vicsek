// compile with: g++ main.cpp Vicsek.cpp nrutil.cpp -o main -O2 -w

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>

#include "nrutil.h"
#include "Vicsek.h"

#define PI 3.1415926

int main (){
    // set parameters
    int n = 10;
    int N = n*n*n;

    double sigma = 1.0;
    double phi = 0.2;

    double D_rot = 0.1;
    double v = 0.03;
    double R = 1.;

    //simulation parameters
    double dt = 1.;// 0.0001; //time-step length
    int Neq   = 2000; // equilibration time
    int Nsim  = 100000; // simulation time
    int Nsave = 1000; // how often printed

    // initialize random number generator
    srand (12345);

    // create the class and the initial configuration
    class Vicsek *system;
    system = new Vicsek(sigma, N, phi, D_rot, v, R);

    // Equilibration
    system->md_equilibration(dt, Neq);

    // Simulation
    system->run_simulation("simulation.out", dt, Nsim, Nsave);

    return 0;
}
