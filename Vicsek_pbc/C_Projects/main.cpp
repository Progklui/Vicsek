// compile with: g++ main.cpp Vicsek.cpp nrutil.cpp -o main -O2 -w

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

#include "ReadData.cpp"

#include "nrutil.h"
#include "Vicsek.h"

#define PI 3.1415926

void getParameterAsArray(double start, double end, double step, int N, double *arr);

int main (int argc, char *argv[]){
    if(argc < 2) {
    	printf("Specify arguments: %s <String>\n", argv[0]);
    	exit(1);
    }

    char *input = argv[1];

    ReadData *newRead = new ReadData(input);

    double *data = newRead->getData();

    // set parameters
    int Ninit = data[0];
    int Nend = data[1];
    int Nincr = data[2];
    double L = data[3];
    bool config = data[17];
    int init = data[18];

    double D_rotInit = data[11];
    double D_rotEnd = data[12];
    double D_rotIncr = data[13];
    int N_D_rot = round((D_rotEnd - D_rotInit) / D_rotIncr);
    double *D_rotArr = new double[N_D_rot];
    getParameterAsArray(D_rotInit, D_rotEnd, D_rotIncr, N_D_rot, D_rotArr);

    double vInit     = data[8];
    double vEnd = data[9];
    double vIncr = data[10];
    int N_v = round((vEnd - vInit) / vIncr);
    double *vArr = new double[N_v];
    getParameterAsArray(vInit, vEnd, vIncr, N_v, vArr);

    double RInit     = data[5];
    double REnd = data[6];
    double RIncr = data[7];
    int N_R = round((REnd - RInit) / RIncr);
    double *RArr = new double[N_R];

    getParameterAsArray(RInit, REnd, RIncr, N_R, RArr);

    //simulation parameters
    double dt = data[4]; // 0.0001; //time-step length
    int Neq   = data[16]; // equilibration time
    int Nsim  = data[14]; // simulation time
    int Nsave = data[15]; // how often printed


    printf("Save Config: %s\n", config ? "true" : "false");

    // initialize random number generator
    srand (12345);

    for(int N = Ninit; N < Nend; N += Nincr) {
    	for(int vN = 0; vN < N_v; vN++) {
    		for(int rN = 0; rN < N_R; rN++) {
    			for(int rD = 0; rD < N_D_rot; rD++) {
			    // settings for directory structure - program automatically creates a sensible directory structure
			    std::string dir_name_phys_params = "N_" + std::to_string(N) + "_L_" + std::to_string(L) + "_v_" + std::to_string(vArr[vN]) + "_R_" + std::to_string(RArr[rN]) + "_D_" + std::to_string(D_rotArr[rD]);
			    std::string dir_name_sim_params  = "Neq_" + std::to_string(Neq) + "_Nsim_" + std::to_string(Nsim) + "_dt_" + std::to_string(dt);
			    std::string dir_name             = "../simulation_results/" + dir_name_phys_params + "/"+ dir_name_sim_params;
			    std::string save_dir_name        = "mkdir -p " + dir_name; // ../simulation_results/test";
			    const int dir_err                = system(save_dir_name.c_str());

			    // create the class and the initial configuration
			    class Vicsek *system;
			    system = new Vicsek(L, N, D_rotArr[rD], vArr[vN], RArr[rN], dt, Nsim, Nsave, &dir_name[0], init);

			    system->md_equilibration(Neq); // Equilibration
			    system->run_simulation(config); // Simulation
			 }
		}
	}
    }

    return 0;
}

void getParameterAsArray(double start, double end, double step, int N, double *arr) {
    //printf("%d\n", N);
    for(int i = 0; i < N; i++) {
    	arr[i] = start + i * step;
    }
}
