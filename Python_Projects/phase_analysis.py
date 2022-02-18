#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 23:12:18 2020

@author: Florian Kluibenschedl
"""
import os, sys

path = os.path.dirname(__file__)
sys.path.append(path)

import calculations as calc # for more involved calculations
import plot
import fit

import numpy as np
import pandas as pd

from uncertainties import unumpy
from uncertainties import ufloat
from uncertainties import ufloat_fromstr

import matplotlib.pyplot as plt

def make_dir_if_not_exists(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

def get_va(folder_list):
    va_mean_array = [] # np.zeros(len(o_folder_path))
    va_std_array  = [] # np.zeros(len(o_folder_path))

    i = 0
    for folder_path in folder_list:
        if os.path.exists(folder_path+'/2_simulation.out'):
            simulation = np.loadtxt(folder_path+'/2_simulation.out', delimiter=' ').T

            va_mean      = np.mean(simulation[1])
            va_std_error = np.std(simulation[1], ddof=1) # / (len(simulation[1])**(0.5))

            va_mean_array.append(va_mean)
            va_std_array.append(va_std_error)

            i += 1
    va_mean_array = np.array(va_mean_array)
    return va_mean_array
# ┌───────────┐
# │ Read data │
# └───────────┘
input_object  = calc.handle_input()
o_folder_path = input_object.get_params()

dim, param_1, param_2 = input_object.get_phase_parameters()

if dim == "1D":
    va_mean_array = get_va(o_folder_path)

    if param_1 == "eta":
        x_data   = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
        x_lab    = r"$D_{rot}$"
        data_lab = r"$\rho = $" + str(input_object.N_1/input_object.L) + r", $R = $" + str(input_object.R_1) + r", $v_a = $" + str(input_object.v_1)
        path     = "../simulation_results/1_phase_transitions/eta_transitions/"
        make_dir_if_not_exists(path)
        img_name = "rho_" + str(input_object.N_1/input_object.L) + "_R_" + str(input_object.R_1) + "_v_a_" + str(input_object.v_1)
    elif param_1 == "rho":
        x_data   = np.arange(input_object.N_1, input_object.N_2, input_object.N_inc)
        x_lab    = r"$N$"
        data_lab = r"$D_{rot} = $" + str(input_object.D_rot_1) + r", $R = $" + str(input_object.R_1) + r", $v_a = $" + str(input_object.v_1)
        path     = "../simulation_results/1_phase_transitions/rho_transition/"
        make_dir_if_not_exists(path)
        img_name = "eta_" + str(input_object.D_rot_1) + "_R_" + str(input_object.R_1) + "_v_a_" + str(input_object.v_1)
    elif param_1 == "r":
        x_data   = np.arange(input_object.R_1, input_object.R_2, input_object.R_inc)
        x_lab    = r"$R$"
        data_lab = r"$\rho = $" + str(input_object.N_1/input_object.L) + r", $\D_{rot} = $" + str(input_object.D_rot_1) + r", $v_a = $" + str(input_object.v_1)
        path     = "../simulation_results/1_phase_transitions/r_transition/"
        make_dir_if_not_exists(path)
        img_name = "rho_" + str(input_object.N_1/input_object.L) + "_eta_" + str(input_object.D_rot_1) + "_v_a_" + str(input_object.v_1)

    phase_1D = plot.without_fit_one_data_line(x=x_data, y=va_mean_array, x_label=x_lab, y_label=r"$v_a$", data_label=data_lab)
    phase_1D.plot_phase(image_name=path + img_name, set_grid=False, set_legend=True)
elif dim == "1D+":
    param_2_values_str = input("Specify values (e. g. 1.34, 2.1): ").split(",")
    # param_2_folders    = [] # np.zeros(len(param_2_values_str))
    x_array     = 0
    va_array    = []
    param_array = []

    for i in range(len(param_2_values_str)):
        if param_2 == "eta":
            D_rot = float(param_2_values_str[i])
            param_array.append(D_rot)
            param_label = r"$D_{rot}$"
            if param_1 == "rho":
                N_array      = np.arange(input_object.N_1, input_object.N_2, input_object.N_inc)
                path_to_file = []
                for N in N_array:
                    path_to_file.append(input_object.get_folder_path(N, input_object.L, input_object.v_1, input_object.R_1, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

                x_array = N_array
                x_label = r"$\rho$"

                path     = "../simulation_results/1_phase_transitions/rho_transition/"
                make_dir_if_not_exists(path)
                img_name = "eta_" + param_2_values_str[0] + "-" + param_2_values_str[len(param_2_values_str)-1] + "_R_" + str(input_object.R_1) + "_v_a_" + str(input_object.v_1)

                plot_label = r"$R = $" + str(input_object.R_1) + r", $v = $" + str(input_object.v_1)

            elif param_1 == "r":
                R_array      = np.arange(input_object.R_1, input_object.R_2, input_object.R_inc)
                path_to_file = []
                for R in R_array:
                    path_to_file.append(input_object.get_folder_path(input_object.N_1, input_object.L, input_object.v_1, R, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

                x_array = R_array
                x_label = r"$R$"

                path     = "../simulation_results/1_phase_transitions/r_transition/"
                make_dir_if_not_exists(path)
                img_name = "rho_" + str(input_object.N_1/input_object.L) + "_eta_" + param_2_values_str[0] + "-" + param_2_values_str[len(param_2_values_str)-1] + "_v_a_" + str(input_object.v_1)

                plot_label = r"$\rho = $" + str(input_object.N_1/input_object.L) + r", $v = $" + str(input_object.v_1)

        elif param_2 == "rho":
            rho = float(param_2_values_str[i])
            param_array.append(rho)
            param_label = r"$\rho$"
            if param_1 == "eta":
                D_rot_array  = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
                path_to_file = []
                for D_rot in D_rot_array:
                    path_to_file.append(input_object.get_folder_path(rho*input_object.L, input_object.L, input_object.v_1, input_object.R_1, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

                x_array = D_rot_array
                x_label = r"$D_{rot}$"

                path     = "../simulation_results/1_phase_transitions/eta_transitions/"
                make_dir_if_not_exists(path)
                img_name = "rho_" + param_2_values_str[0] + "-" + param_2_values_str[len(param_2_values_str)-1] + "_R_" + str(input_object.R_1) + "_v_a_" + str(input_object.v_1)

                plot_label = r"$R = $" + str(input_object.R_1) + r", $v = $" + str(input_object.v_1)

            elif param_1 == "r":
                R_array      = np.arange(input_object.R_1, input_object.R_2, input_object.R_inc)
                path_to_file = []
                for R in R_array:
                    path_to_file.append(input_object.get_folder_path(rho*input_object.L, input_object.L, input_object.v_1, R, input_object.D_rot_1, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

                x_array = R_array
                x_label = r"$R$"

                path     = "../simulation_results/1_phase_transitions/r_transition/"
                make_dir_if_not_exists(path)
                img_name = "rho_" + param_2_values_str[0] + "-" + param_2_values_str[len(param_2_values_str)-1] + "_eta_" + str(input_object.D_rot_1) + "_v_a_" + str(input_object.v_1)

                plot_label = r"$D_{rot} = $" + str(input_object.D_rot_1) + r", $v = $" + str(input_object.v_1)

        elif param_2 == "r":
            R = float(param_2_values_str[i])
            param_array.append(D_rot)
            param_label = r"$R$"
            if param_1 == "rho":
                N_array     = np.arange(input_object.N_1, input_object.N_2, input_object.N_inc)
                path_to_file = []
                for N in N_array:
                    path_to_file.append(input_object.get_folder_path(N, input_object.L, input_object.v_1, R, input_object.D_rot_1, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

                x_array = N_array
                x_label = r"$\rho$"

                path     = "../simulation_results/1_phase_transitions/rho_transition/"
                make_dir_if_not_exists(path)
                img_name = "eta_" + str(input_object.D_rot_1) + "_R_" + param_2_values_str[0] + "-" + param_2_values_str[len(param_2_values_str)-1] + "_v_a_" + str(input_object.v_1)

                plot_label = r"$D_{rot} = $" + str(input_object.D_rot_1) + r", $v = $" + str(input_object.v_1)

            elif param_1 == "eta":
                D_rot_array  = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
                path_to_file = []
                for D_rot in D_rot_array:
                    path_to_file.append(input_object.get_folder_path(input_object.N_1, input_object.L, input_object.v_1, R, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

                x_array = D_rot_array
                x_label = r"$D_{rot}$"

                path     = "../simulation_results/1_phase_transitions/eta_transitions/"
                make_dir_if_not_exists(path)
                img_name = "rho_" + str(input_object.N_1/input_object.L) + "_R_" + param_2_values_str[0] + "-" + param_2_values_str[len(param_2_values_str)-1] + "_v_a_" + str(input_object.v_1)
                plot_label = r"$\rho = $" + str(input_object.N_1/input_object.L) + r", $v_a = $" + str(input_object.v_1)

    va_array    = np.array(va_array)
    param_array = np.array(param_array)

    plot_object = plot.phase_analysis(x=x_array, y=va_array, params=param_array, x_label=x_label, y_label=r"$v_a$", plot_label=plot_label)
    plot_object.plot_phases(param_label=param_label, image_name=path + img_name, set_grid=False, set_legend=True)

elif dim == "2D":
    va_array = []
    if param_1 == "rho":
        x_label = r"$\rho$"
        N_array      = np.arange(input_object.N_1, input_object.N_2, input_object.N_inc)
        x_data = N_array
        for N in N_array:
            if param_2 == "eta":
                y_label = r"$D_{rot}$"
                plot_label = r"$R = $" + str(input_object.R_1) + r", $v_a = $" + str(input_object.v_1)

                path     = "../simulation_results/1_phase_transitions/rho_transition/"
                make_dir_if_not_exists(path)
                img_name = "R_" + str(input_object.R_1) + "_v_a_" + str(input_object.v_1)

                D_rot_array  = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
                y_data = D_rot_array
                path_to_file = []
                for D_rot in D_rot_array:
                    path_to_file.append(input_object.get_folder_path(N, input_object.L, input_object.v_1, input_object.R_1, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)
            elif param_2 == "r":
                y_label = r"$R$"
                plot_label = r"$D_{rot} = $" + str(input_object.D_rot_1) + r", $v_a = $" + str(input_object.v_1)

                path     = "../simulation_results/1_phase_transitions/rho_transition/"
                make_dir_if_not_exists(path)
                img_name = "eta_" + str(input_object.D_rot_1) + "_v_a_" + str(input_object.v_1)

                R_array      = np.arange(input_object.R_1, input_object.R_2, input_object.R_inc)
                y_data = R_array
                path_to_file = []
                for R in R_array:
                    path_to_file.append(input_object.get_folder_path(N, input_object.L, input_object.v_1, R, input_object.D_rot_1, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

    elif param_1 == "eta":
        x_label = r"$D_{rot}$"
        D_rot_array  = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
        x_data = D_rot_array
        for D_rot in D_rot_array:
            if param_2 == "rho":
                y_label = r"$\rho$"
                plot_label = r"$R = $" + str(input_object.R_1) + r", $v_a = $" + str(input_object.v_1)

                path     = "../simulation_results/1_phase_transitions/eta_transitions/"
                make_dir_if_not_exists(path)
                img_name = "R_" + str(input_object.R_1) + "_v_a_" + str(input_object.v_1)

                N_array      = np.arange(input_object.N_1, input_object.N_2, input_object.N_inc)
                y_data = N_array
                path_to_file = []
                for N in N_array:
                    path_to_file.append(input_object.get_folder_path(N, input_object.L, input_object.v_1, input_object.R_1, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)
            elif param_2 == "r":
                y_label = r"$R$"
                plot_label = r"$\rho = $" + str(input_object.N_1/input_object.L) + r", $v_a = $" + str(input_object.v_1)

                path     = "../simulation_results/1_phase_transitions/eta_transitions/"
                make_dir_if_not_exists(path)
                img_name = "rho_" + str(input_object.N_1/input_object.L) + "_v_a_" + str(input_object.v_1)

                R_array      = np.arange(input_object.R_1, input_object.R_2, input_object.R_inc)
                y_data = R_array
                path_to_file = []
                for R in R_array:
                    path_to_file.append(input_object.get_folder_path(input_object.N_1, input_object.L, input_object.v_1, R, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

    elif param_1 == "r":
        x_label = r"$R$"
        R_array      = np.arange(input_object.R_1, input_object.R_2, input_object.R_inc)
        x_data = R_array
        for R in R_array:
            if param_2 == "rho":
                y_label = r"$\rho$"
                plot_label = r"$D_{rot} = $" + str(input_object.D_rot_1) + r", $v_a = $" + str(input_object.v_1)

                path     = "../simulation_results/1_phase_transitions/r_transition/"
                make_dir_if_not_exists(path)
                img_name = "eta_" + str(input_object.D_rot_1) + "_v_a_" + str(input_object.v_1)

                N_array     = np.arange(input_object.N_1, input_object.N_2, input_object.N_inc)
                y_data = N_array
                path_to_file = []
                for N in N_array:
                    path_to_file.append(input_object.get_folder_path(N, input_object.L, input_object.v_1, R, input_object.D_rot_1, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

            elif param_2 == "eta":
                y_label = r"$D_{rot}$"
                plot_label = r"$\rho = $" + str(input_object.N_1/input_object.L) + r", $v_a = $" + str(input_object.v_1)

                path     = "../simulation_results/1_phase_transitions/r_transition/"
                make_dir_if_not_exists(path)
                img_name = "rho_" + str(input_object.N_1/input_object.L) + "_v_a_" + str(input_object.v_1)

                D_rot_array  = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
                y_data = D_rot_array
                path_to_file = []
                for D_rot in D_rot_array:
                    path_to_file.append(input_object.get_folder_path(input_object.N_1, input_object.L, input_object.v_1, R, D_rot, input_object.Neq, input_object.Nsim, input_object.dt))
                va_mean_array = get_va(path_to_file)
                va_array.append(va_mean_array)

    print(x_data)
    print(y_data)

    X, Y = np.meshgrid(x_data, y_data)
    ax = plt.subplot(111)
    ax.contourf(X, Y, np.array(va_array).T)
    plt.show()
