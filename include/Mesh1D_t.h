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
};

#endif