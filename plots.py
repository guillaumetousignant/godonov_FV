# MCG 5136 Assignment 3
# Guillaume Tousignant, 0300151859
# November 27th, 2020

import matplotlib.pyplot as plt
import numpy as np
import os
import re
from pathlib import Path

def find_problem_files(filter):
    times = []
    problem_numbers = []
    x_arrays = []
    ux_arrays = []
    rho_arrays = []
    p_arrays = []
    mach_arrays = []
    T_arrays = []

    t_finder = re.compile(r"SOLUTIONTIME = [-+]?\d*\.?\d+")
    I_finder = re.compile(r"I= \d*")
    N_finder = re.compile(r"N= \d*")
    problem_finder = re.compile(r"Problem \d*")

    # Input from all the output_X.dat files
    filenames = [f for f in os.listdir(os.path.join(os.getcwd(), 'data')) if os.path.isfile(os.path.join(os.getcwd(), 'data', f)) and "output_" in f and filter in f and f.endswith(".dat")]
    for filename in filenames:
        with open(os.path.join(os.getcwd(), 'data', filename), 'r') as file:
            lines = file.readlines()
            t_match = t_finder.search(lines[2])
            times.append(float(t_match.group(0)[15:]))
            N_points_match = I_finder.search(lines[2])
            N_elements_match = N_finder.search(lines[0])
            N_points = int(N_points_match.group(0)[3:])
            N_elements = int(N_elements_match.group(0)[3:])
            problem_match = problem_finder.search(lines[0])
            problem_numbers.append(int(problem_match.group(0)[8:]))

            x_arrays.append(np.zeros(N_points))
            ux_arrays.append(np.zeros(N_points))
            rho_arrays.append(np.zeros(N_points))
            p_arrays.append(np.zeros(N_points))
            mach_arrays.append(np.zeros(N_points))
            T_arrays.append(np.zeros(N_points))

            for i in range(N_points):
                numbers = lines[i+3].split()
                x_arrays[-1][i] = float(numbers[0])
                ux_arrays[-1][i] = float(numbers[1])
                rho_arrays[-1][i] = float(numbers[2])
                p_arrays[-1][i] = float(numbers[3])
                mach_arrays[-1][i] = float(numbers[4])
                T_arrays[-1][i] = float(numbers[5])

    return zip(*sorted(zip(problem_numbers, filenames, times, x_arrays, ux_arrays, rho_arrays, p_arrays, mach_arrays, T_arrays)))

problem_numbers, filenames, times, x_arrays, ux_arrays, rho_arrays, p_arrays, mach_arrays, T_arrays = find_problem_files("_exact")

# Plotting
save_path = Path.cwd() / "figures"
save_path.mkdir(parents=True, exist_ok=True)

for i in range(len(filenames)):
    u_fig, u_ax = plt.subplots(1, 1)
    u_ax.plot(x_arrays[i], ux_arrays[i])

    u_ax.grid()
    u_ax.set_ylabel('$U_x$ [$m/s$]')
    u_ax.set_xlabel('x [m]')
    u_ax.set_title(f"Problem {problem_numbers[i]}, $U_x$ along x at t = {times[i]}s")
    u_fig.canvas.set_window_title(f'Problem {problem_numbers[i]} $U_x$')
    u_fig.savefig(save_path / f"problem_{problem_numbers[i]}_u.svg", format='svg', transparent=True)
    u_fig.savefig(save_path / f"problem_{problem_numbers[i]}_u.png", format='png', transparent=True)

    rho_fig, rho_ax = plt.subplots(1, 1)
    rho_ax.plot(x_arrays[i], rho_arrays[i])

    rho_ax.grid()
    rho_ax.set_ylabel(r'$\rho$ [$kg/m^3$]')
    rho_ax.set_xlabel('x [m]')
    rho_ax.set_title(f"Problem {problem_numbers[i]}, $\\rho$ along x at t = {times[i]}s")
    rho_fig.canvas.set_window_title(f'Problem {problem_numbers[i]} $\\rho$')
    rho_fig.savefig(save_path / f"problem_{problem_numbers[i]}_rho.svg", format='svg', transparent=True)
    rho_fig.savefig(save_path / f"problem_{problem_numbers[i]}_rho.png", format='png', transparent=True)

    p_fig, p_ax = plt.subplots(1, 1)
    p_ax.plot(x_arrays[i], p_arrays[i]/1000)

    p_ax.grid()
    p_ax.set_ylabel('p [kPa]')
    p_ax.set_xlabel('x [m]')
    p_ax.set_title(f"Problem {problem_numbers[i]}, p along x at t = {times[i]}s")
    p_fig.canvas.set_window_title(f'Problem {problem_numbers[i]} p')
    p_fig.savefig(save_path / f"problem_{problem_numbers[i]}_p.svg", format='svg', transparent=True)
    p_fig.savefig(save_path / f"problem_{problem_numbers[i]}_p.png", format='png', transparent=True)

    mach_fig, mach_ax = plt.subplots(1, 1)
    mach_ax.plot(x_arrays[i], mach_arrays[i])

    mach_ax.grid()
    mach_ax.set_ylabel('M [-]')
    mach_ax.set_xlabel('x [m]')
    mach_ax.set_title(f"Problem {problem_numbers[i]}, M along x at t = {times[i]}s")
    mach_fig.canvas.set_window_title(f'Problem {problem_numbers[i]} M')
    mach_fig.savefig(save_path / f"problem_{problem_numbers[i]}_mach.svg", format='svg', transparent=True)
    mach_fig.savefig(save_path / f"problem_{problem_numbers[i]}_mach.png", format='png', transparent=True)

    T_fig, T_ax = plt.subplots(1, 1)
    T_ax.plot(x_arrays[i], T_arrays[i])

    T_ax.grid()
    T_ax.set_ylabel('T [K]')
    T_ax.set_xlabel('x [m]')
    T_ax.set_title(f"Problem {problem_numbers[i]}, T along x at t = {times[i]}s")
    T_fig.canvas.set_window_title(f'Problem {problem_numbers[i]} T')
    T_fig.savefig(save_path / f"problem_{problem_numbers[i]}_T.svg", format='svg', transparent=True)
    T_fig.savefig(save_path / f"problem_{problem_numbers[i]}_T.png", format='png', transparent=True)

plt.show()