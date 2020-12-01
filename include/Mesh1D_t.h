#ifndef MESH1D_T_H
#define MESH1D_T_H

class Mesh1D_t { 
    public: 
        Mesh1D_t(int n_cells);
        ~Mesh1D_t();

        int n_cells_;
        double* rho_;
        double* u_;
        double* p_;
        double* x_;

        void initial_conditions(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double discontinuity);
};

#endif