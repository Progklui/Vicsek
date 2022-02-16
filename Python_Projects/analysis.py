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

# ┌───────────┐
# │ Read data │
# └───────────┘
input_object = calc.handle_input()
folder_path  = input_object.get_params()

equilibration = np.loadtxt(folder_path+'/1_equilibration.out', delimiter=' ').T # , skiprows=0
simulation    = np.loadtxt(folder_path+'/2_simulation.out', delimiter=' ').T

# ┌──────────────────────────┐
# │ Calculation & or Fitting │
# └──────────────────────────┘



# ┌───────────┐
# │ Plot data │
# └───────────┘

print_traj = input_object.get_trajectory()
if print_traj == True:
    equi_trajectory = plot.without_fit_one_data_line(x=equilibration[0], y=equilibration[1],
                                                     x_label="time", y_label=r"$v_a$", data_label="equilibration")
    equi_trajectory.scatter(image_name=folder_path + "/1_equilibration_trajectory", set_grid=True, set_legend=True)

    sim_trajectory = plot.without_fit_one_data_line(x=simulation[0], y=simulation[1],
                                                    x_label="time", y_label=r"$v_a$", data_label="simulation")
    sim_trajectory.scatter(image_name=folder_path + "/2_simulation_trajectory", set_grid=True, set_legend=True)

cont = True
while cont == True:
    config_file_name, cont = input_object.get_configuration()

    x     = np.array(pd.read_csv(folder_path + config_file_name, usecols=[0], delimiter=" "))[:, -1]
    y     = np.array(pd.read_csv(folder_path + config_file_name, usecols=[1], delimiter=" "))[:, -1]
    theta = np.array(pd.read_csv(folder_path + config_file_name, usecols=[2], delimiter=" "))[:, -1]

    config_object = plot.configuration(x=x, y=y, theta=theta, x_label=r"x", y_label=r"y", data_label=r"test")
    config_object.plot(image_name=folder_path + "/3_config", set_grid=False, set_legend=True)
# ┌─────────────┐
# │ Export data │
# └─────────────┘
