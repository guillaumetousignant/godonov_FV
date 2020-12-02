#include "Solver_t.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <iomanip>

namespace fs = std::filesystem;

Solver_t::Solver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double end_time, double discontinuity, int n_points, int problem_number) :
        gamma_{1.4, 1.4},
        R_(8.31446261815324),
        a_{std::sqrt(gamma_[0] * p_L/rho_L), std::sqrt(gamma_[1] * p_R/rho_R)},
        u_{u_L, u_R},
        p_{p_L, p_R},
        end_time_(end_time),
        discontinuity_(discontinuity),
        n_points_(n_points),
        problem_number_(problem_number),
        x_{x_L, x_R} {


}

Solver_t::~Solver_t() {}

void Solver_t::write_file_data(int N_points, double time, const std::vector<double> &rho, const std::vector<double> &u, const std::vector<double> &p, const std::vector<double> &x, const std::vector<double> &mach, const std::vector<double> &T, int problem_number, std::string suffix) {
    std::stringstream ss;
    std::ofstream file;

    fs::path save_dir = fs::current_path() / "data";
    fs::create_directory(save_dir);

    ss << "output_" << problem_number << suffix << ".dat";
    file.open(save_dir / ss.str());

    file << "TITLE = \"Problem " << problem_number << " at t= " << time << "\"" << std::endl;
    file << "VARIABLES = \"X\", \"U_x\", \"rho\", \"p\", \"mach\", \"T\"" << std::endl;
    file << "ZONE T= \"Zone     1\",  I= " << N_points << ",  J= 1,  DATAPACKING = POINT, SOLUTIONTIME = " << time << std::endl;

    for (int i = 0; i < N_points; ++i) {
        file << std::setw(12) << x[i] << " " << std::setw(12) << u[i] << " " << std::setw(12) << rho[i] << " " << std::setw(12) << p[i] << " " << std::setw(12) << mach[i] << " " << std::setw(12) << T[i] << std::endl;
    }

    file.close();
}