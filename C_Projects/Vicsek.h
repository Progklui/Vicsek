class Vicsek {
    public:
	      Vicsek(double, int, double, double, double, double);

	      double get_size();
        double get_packing_fraction();
        double calculate_va();

        void md_step_vicsek(double);
        void md_equilibration(char * str, double, int);
        void run_simulation(char * str, double, int, int);

    private:
        int N;
        int n;
        double L;
        double sigma;
        double sigma2;
        double D_rot; // rotational diffusion coefficient
        double v; // magnitude of the velocity of each particle
        double R; // radius that is taken into account to calculate angle

        // positions x,y,z
        double * x;
        double * y;

        double * x_init;
        double * y_init;

        //angle vector accounting for orientation
        double * theta;

        // private functions
        void pbc(double &x, double &y);
        double calculate_mean_angle(int);
        void store_configuration(double);
};
