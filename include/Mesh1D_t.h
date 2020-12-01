#ifndef MESH1D_T_H
#define MESH1D_T_H

class Mesh1D_t { 
    public: 
        Mesh1D_t(int n_cells);
        ~Mesh1D_t();

        int n_cells_;
        double* a_;
        double* u_;
        double* p_;
        double* x_;
        double* gamma_;
        double delta_x_;
        int n_faces_;
        double* F_1_; // ρu
        double* F_2_; // ρu² + p
        double* F_3_; // u(γp/(γ-1) + ρu²/2)

        void initial_conditions(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double gamma_L, double gamma_R, double discontinuity);
};

#endif