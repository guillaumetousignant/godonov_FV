#ifndef MESH1D_T_H
#define MESH1D_T_H

#include <vector>

class Mesh1D_t { 
    public: 
        Mesh1D_t(int n_cells, double delta_x);
        ~Mesh1D_t();

        int n_cells_;
        std::vector<double> a_;
        std::vector<double> u_;
        std::vector<double> p_;
        std::vector<double> x_;
        std::vector<double> gamma_;
        double delta_x_;
        int n_faces_;
        std::vector<double> F_1_; // ρu
        std::vector<double> F_2_; // ρu² + p
        std::vector<double> F_3_; // u(γp/(γ-1) + ρu²/2)

        void initial_conditions(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double gamma_L, double gamma_R, double discontinuity);
};

#endif