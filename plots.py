# MCG 5136 Assignment 3
# Guillaume Tousignant, 0300151859
# November 27th, 2020

import matplotlib.pyplot as plt
import numpy as np
import os
import re
from pathlib import Path

class Result:
    pass

def find_problem_files(filter):
    result = Result()
    result.times = []
    result.problem_numbers = []
    result.x_arrays = []
    result.ux_arrays = []
    result.rho_arrays = []
    result.p_arrays = []
    result.mach_arrays = []
    result.T_arrays = []
    result.N_arrays = []

    t_finder = re.compile(r"SOLUTIONTIME = [-+]?\d*\.?\d+")
    I_finder = re.compile(r"I= \d*")
    N_finder = re.compile(r"N= \d*")
    problem_finder = re.compile(r"Problem \d*")

    # Input from all the output_X.dat files
    result.filenames = [f for f in os.listdir(os.path.join(os.getcwd(), 'data')) if os.path.isfile(os.path.join(os.getcwd(), 'data', f)) \
                            and "output_" in f \
                            and filter in f \
                            and f.endswith(".dat")]

    for filename in result.filenames:
        with open(os.path.join(os.getcwd(), 'data', filename), 'r') as file:
            lines = file.readlines()
            t_match = t_finder.search(lines[2])
            result.times.append(float(t_match.group(0)[15:]))
            N_points_match = I_finder.search(lines[2])
            N_elements_match = N_finder.search(lines[0])
            N_points = int(N_points_match.group(0)[3:])
            problem_match = problem_finder.search(lines[0])
            result.problem_numbers.append(int(problem_match.group(0)[8:]))
            result.N_arrays.append(int(N_elements_match.group(0)[3:]))

            result.x_arrays.append(np.zeros(N_points))
            result.ux_arrays.append(np.zeros(N_points))
            result.rho_arrays.append(np.zeros(N_points))
            result.p_arrays.append(np.zeros(N_points))
            result.mach_arrays.append(np.zeros(N_points))
            result.T_arrays.append(np.zeros(N_points))

            for i in range(N_points):
                numbers = lines[i+3].split()
                result.x_arrays[-1][i] = float(numbers[0])
                result.ux_arrays[-1][i] = float(numbers[1])
                result.rho_arrays[-1][i] = float(numbers[2])
                result.p_arrays[-1][i] = float(numbers[3])
                result.mach_arrays[-1][i] = float(numbers[4])
                result.T_arrays[-1][i] = float(numbers[5])

    # Put in problem order
    result.problem_numbers, result.filenames, result.times, result.x_arrays, result.ux_arrays, result.rho_arrays, result.p_arrays, result.mach_arrays, result.T_arrays, result.N_arrays = \
        zip(*sorted(zip(result.problem_numbers, result.filenames, result.times, result.x_arrays, result.ux_arrays, result.rho_arrays, result.p_arrays, result.mach_arrays, result.T_arrays, result.N_arrays)))

    return result

def find_problem_files_N(filter, N_vec):
    result = Result()
    result.problem_numbers = []
    result.filenames = []
    result.times = []
    result.x_arrays = []
    result.ux_arrays = []
    result.rho_arrays = []
    result.p_arrays = []
    result.mach_arrays = []
    result.T_arrays = []
    result.N_arrays = []
    result.N = N_vec

    for i in result.N:
        result_N = find_problem_files(f"{filter}_N{i}.")
        result.problem_numbers.append(result_N.problem_numbers)
        result.filenames.append(result_N.filenames)
        result.times.append(result_N.times)
        result.x_arrays.append(result_N.x_arrays)
        result.ux_arrays.append(result_N.ux_arrays)
        result.rho_arrays.append(result_N.rho_arrays)
        result.p_arrays.append(result_N.p_arrays)
        result.mach_arrays.append(result_N.mach_arrays)
        result.T_arrays.append(result_N.T_arrays)
        result.N_arrays.append(result_N.N_arrays)

    return result

def generate_figures(save_path, result_exact, result, title, suffix):
    for i in range(len(result_exact.problem_numbers)):
        u_fig, u_ax = plt.subplots(1, 1)
        u_ax.plot(result_exact.x_arrays[i], result_exact.ux_arrays[i], label="Exact solution")
        for j in range(len(result.N)):
            u_ax.plot(result.x_arrays[j][i], result.ux_arrays[j][i], label=f"{title} fluxes, N = {N[j]}")

        u_ax.grid()
        u_ax.set_ylabel('$U_x$ [$m/s$]')
        u_ax.set_xlabel('x [m]')
        u_ax.set_title(f"Problem {result_exact.problem_numbers[i]}, {title} fluxes, $U_x$ along x at t = {result_exact.times[i]}s")
        u_ax.legend()
        u_fig.canvas.set_window_title(f'Problem {result_exact.problem_numbers[i]}, {title} fluxes U_x')
        u_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_u_{suffix}.svg", format='svg', transparent=True)
        u_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_u_{suffix}.png", format='png', transparent=True)

        rho_fig, rho_ax = plt.subplots(1, 1)
        rho_ax.plot(result_exact.x_arrays[i], result_exact.rho_arrays[i], label="Exact solution")
        for j in range(len(result.N)):
            rho_ax.plot(result.x_arrays[j][i], result.rho_arrays[j][i], label=f"{title} fluxes, N = {N[j]}")

        rho_ax.grid()
        rho_ax.set_ylabel(r'$\rho$ [$kg/m^3$]')
        rho_ax.set_xlabel('x [m]')
        rho_ax.set_title(f"Problem {result_exact.problem_numbers[i]}, {title} fluxes, $\\rho$ along x at t = {result_exact.times[i]}s")
        rho_ax.legend()
        rho_fig.canvas.set_window_title(f'Problem {result_exact.problem_numbers[i]}, {title} fluxes rho')
        rho_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_rho_{suffix}.svg", format='svg', transparent=True)
        rho_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_rho_{suffix}.png", format='png', transparent=True)

        p_fig, p_ax = plt.subplots(1, 1)
        p_ax.plot(result_exact.x_arrays[i], result_exact.p_arrays[i]/1000, label="Exact solution")
        for j in range(len(result.N)):
            p_ax.plot(result.x_arrays[j][i], result.p_arrays[j][i]/1000, label=f"{title} fluxes, N = {N[j]}")

        p_ax.grid()
        p_ax.set_ylabel('p [kPa]')
        p_ax.set_xlabel('x [m]')
        p_ax.set_title(f"Problem {result_exact.problem_numbers[i]}, {title} fluxes, p along x at t = {result_exact.times[i]}s")
        p_ax.legend()
        p_fig.canvas.set_window_title(f'Problem {result_exact.problem_numbers[i]}, {title} fluxes p')
        p_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_p_{suffix}.svg", format='svg', transparent=True)
        p_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_p_{suffix}.png", format='png', transparent=True)

        mach_fig, mach_ax = plt.subplots(1, 1)
        mach_ax.plot(result_exact.x_arrays[i], result_exact.mach_arrays[i], label="Exact solution")
        for j in range(len(result.N)):
            mach_ax.plot(result.x_arrays[j][i], result.mach_arrays[j][i], label=f"{title} fluxes, N = {N[j]}")

        mach_ax.grid()
        mach_ax.set_ylabel('M [-]')
        mach_ax.set_xlabel('x [m]')
        mach_ax.set_title(f"Problem {result_exact.problem_numbers[i]}, {title} fluxes, M along x at t = {result_exact.times[i]}s")
        mach_ax.legend()
        mach_fig.canvas.set_window_title(f'Problem {result_exact.problem_numbers[i]}, {title} fluxes M')
        mach_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_mach_{suffix}.svg", format='svg', transparent=True)
        mach_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_mach_{suffix}.png", format='png', transparent=True)

        T_fig, T_ax = plt.subplots(1, 1)
        T_ax.plot(result_exact.x_arrays[i], result_exact.T_arrays[i], label="Exact solution")
        for j in range(len(result.N)):
            T_ax.plot(result.x_arrays[j][i], result.T_arrays[j][i], label=f"{title} fluxes, N = {N[j]}")

        T_ax.grid()
        T_ax.set_ylabel('T [K]')
        T_ax.set_xlabel('x [m]')
        T_ax.set_title(f"Problem {result_exact.problem_numbers[i]}, {title} fluxes, T along x at t = {result_exact.times[i]}s")
        T_ax.legend()
        T_fig.canvas.set_window_title(f'Problem {result_exact.problem_numbers[i]}, {title} fluxes T')
        T_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_T_{suffix}.svg", format='svg', transparent=True)
        T_fig.savefig(save_path / f"problem_{result_exact.problem_numbers[i]}_T_{suffix}.png", format='png', transparent=True)


N = [100, 200, 500, 1000]

result_exact = find_problem_files("exact")
result_riemann = find_problem_files_N("riemann", N)
#result_roe = find_problem_files_N("roe", N)
#result_roe_entropy = find_problem_files_N("roe_entropy", N)
#result_hlle = find_problem_files_N("hlle", N)


# Plotting
save_path = Path.cwd() / "figures"
save_path.mkdir(parents=True, exist_ok=True)

# Plotting everything at the same time will stack overflow.
generate_figures(save_path, result_exact, result_riemann, "Riemann problem", "riemann")
#generate_figures(save_path, result_exact, result_roe, "Roe", "roe")
#generate_figures(save_path, result_exact, result_roe_entropy, "Roe entropy fix", "roe_entropy")
#generate_figures(save_path, result_exact, result_hlle, "HLLE", "hlle")

plt.show()