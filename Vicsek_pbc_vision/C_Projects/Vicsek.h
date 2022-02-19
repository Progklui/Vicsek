class Vicsek {
    public:
	      Vicsek(double, int, double, double, double, double, int, int, char *, int, double);

        double calculate_va();

        void md_step_vicsek(double);
        void md_equilibration(int);
        void run_simulation(bool);

    private:
        double dt;
        int Nsim;
        int Nsave;

        int N;
        int n;
        double L;
        double D_rot; // rotational diffusion coefficient
        double v; // magnitude of the velocity of each particle
        double R; // radius that is taken into account to calculate angle

        double phi_vision;

        // positions x,y,z
        double * x;
        double * y;

        double * x_init;
        double * y_init;

        //angle vector accounting for orientation
        double * theta;

        char * dir_name;

        // private functions
        void pbc(double &x, double &y);
        void rbc(double &x, double &y);
        double calculate_mean_angle(int);
        void store_configuration(double);
        void save_simulation_params();
        void init_uniform();
        void init_random();
        void check_vision(int, int);
};
