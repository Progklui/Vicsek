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

# ┌───────────┐
# │ Read data │
# └───────────┘
input_object   = calc.handle_input()
o_folder_path  = input_object.get_params()

other_va = []
critical_exp_va_mean = [] #np.zeros(int(len(o_folder_path)/2))

eta_array_all = np.arange(input_object.D_rot_1, input_object.D_rot_2, input_object.D_rot_inc)
eta_array_all = (2*eta_array_all)**0.5
eta_array     = []
other_eta = []

critical_eta = float(input("Input critical eta: "))
va_threshold = float(input("va threshold: "))

i = 0
for folder_path in o_folder_path:
    equilibration = np.loadtxt(folder_path+'/1_equilibration.out', delimiter=' ').T # , skiprows=0
    simulation    = np.loadtxt(folder_path+'/2_simulation.out', delimiter=' ').T

    va_mean = np.mean(simulation[1])
    if eta_array_all[i] < critical_eta:
        if va_mean > va_threshold:# and va_mean > 0.2:
            critical_exp_va_mean.append(va_mean)
            eta_array.append(eta_array_all[i])
        else:
            other_va.append(va_mean)
            other_eta.append(eta_array_all[i])

    i += 1
    if i == 100:
        break

other_va  = np.array(other_va)
other_eta = np.array(other_eta)

critical_exp_va_mean = np.array(critical_exp_va_mean)
eta_array            = np.array(eta_array)

x_data = (critical_eta - eta_array)/critical_eta
y_data = critical_exp_va_mean

print(x_data)
print(y_data)

range_low = 0

func_model = fit.model
fit_object = fit.fit_two_params(function=func_model.linear, x=np.log(x_data), y=y_data, y_error=y_data)
fit_data, *params = fit_object.unweighted(start=(1,2))

plt.scatter(x_data, y_data, marker="x", color="black", label=r"data points for fit")
plt.plot(np.exp(fit_data[0]), fit_data[1], color="red", label=r"linear fit")

plt.scatter((critical_eta - other_eta)/critical_eta, other_va, marker="x", label=r"not fitted")

plt.axhline(y=1, color="grey", linestyle="dashed", linewidth=1)

plt.xlabel(r"$(\sqrt{2D_{rot}^{c}\Delta t} - \sqrt{2D_{rot}\Delta t})/\sqrt{2D_{rot}^{c}\Delta t}$")
plt.ylabel(r"$v_a$")

plt.legend(loc="upper left")
plt.xscale('log')

plt.savefig("critical_exponents.png", dpi=800)
plt.show()
