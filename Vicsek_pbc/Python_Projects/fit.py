#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 23:12:18 2020

@author: Florian Kluibenschedl
"""
import os, sys

path = os.path.dirname(__file__)
sys.path.append(path)

import numpy as np

from scipy.optimize import curve_fit

from uncertainties import unumpy
from uncertainties import ufloat

# ┌───────────────┐
# │ Define models │
# └───────────────┘
class model:    
    def linear(x, a, b):
        return a * x + b

    def quadratic(x, a, b, c):
        return a * x**2 + b * x + c

# ┌─────────────────────────┐
# │ Fit routines - 2 params │
# └─────────────────────────┘
class fit_two_params:
    def __init__(self, function, x, y, y_error):
        self.function       = function
        self.x              = x
        self.y              = y
        self.y_error        = y_error
        self.y_error_scaled = 0
        self.red_chi_squ    = 1
        self.r_squ          = 1
        self.dof            = 0
        
    def get_unweighted_model(self, start):
        p_opt, p_cov = curve_fit(self.function, xdata=self.x, ydata=self.y, 
                                 p0=start, absolute_sigma=False, maxfev=10000)
    
        p_err = np.sqrt(np.diag(p_cov))

        fit_a = ufloat(p_opt[0], p_err[0])
        fit_b = ufloat(p_opt[1], p_err[1])
    
        ## Calculation of R^2
        residuals = self.y - self.function(self.x, *p_opt) 
        ss_res = np.sum(residuals**2) 
        ss_tot = np.sum((self.y - np.mean(self.y))**2) 

        self.r_squ = 1 - (ss_res / ss_tot)
        self.dof   = len(self.x) - 2
    
        errors_scaled = np.sum(((self.y - self.function(self.x, *p_opt)))**2)**0.5
        self.y_error_scaled = unumpy.std_devs(unumpy.uarray(unumpy.nominal_values(self.y), errors_scaled))
        
        print("Unweighted fit with two parameters:")
        print("1st parameter =", fit_a, " unit")
        print("2nd parameter =", fit_b, " unit")
        print("R^2 (2 parameter model) =", round(self.r_squ, 4), ", dof =", self.dof)
        print(" ")
    
        return fit_a, fit_b

    def unweighted(self, start):
        fit_a, fit_b = self.get_unweighted_model(start)
    
        x_fit = np.linspace(self.x[0], self.x[len(self.x)-1], 1000)
        y_fit = self.function(x_fit, fit_a.n, fit_b.n)
        
        return np.array([x_fit, y_fit]), fit_a, fit_b

    def get_weighted_model(self, start):
        p_opt, p_cov = curve_fit(self.function, xdata=self.x, ydata=self.y, 
                                 sigma=self.y_error, p0=start, absolute_sigma=False, maxfev=10000)
    
        p_err = np.sqrt(np.diag(p_cov))

        fit_a = ufloat(p_opt[0], p_err[0])
        fit_b = ufloat(p_opt[1], p_err[1])
        
        ## Calculation of R^2
        residuals = self.y - self.function(self.x, *p_opt) 
        ss_res = np.sum(residuals**2) 
        ss_tot = np.sum((self.y - np.mean(self.y))**2) 
        
        self.r_squ = 1 - (ss_res / ss_tot)
        self.dof = len(self.x) - 2
        
        self.red_chi_squ = np.sum(((self.y - self.function(self.x, *p_opt)) / self.y_error)**2) / self.dof
        
        errors_scaled = np.sum(((self.y - self.function(self.x, *p_opt)))**2)**0.5
        self.y_error_scaled = unumpy.std_devs(unumpy.uarray(unumpy.nominal_values(self.y), errors_scaled))
        
        print("Weighted fit with two parameters:")
        print("1st parameter =", fit_a, " unit")
        print("2nd parameter =", fit_b, " unit")
        print("R^2 (2 parameter model) =", round(self.r_squ,4), ", dof =", self.dof)
        print("red Chi^2 (2 parameter model) =", round(self.red_chi_squ,4), ", dof =", self.dof)
        print(" ")
    
        return fit_a, fit_b

    def weighted(self, start):
        fit_a, fit_b = self.get_weighted_model(start)
    
        x_fit = np.linspace(self.x[0], self.x[len(self.x)-1], 1000)
        y_fit = self.function(x_fit, fit_a.n, fit_b.n)
        
        return np.array([x_fit, y_fit]), fit_a, fit_b


# ┌─────────────────────────┐
# │ Fit routines - 3 params │
# └─────────────────────────┘
class fit_three_params:
    def __init__(self, function, x, y, y_error):
        self.function       = function
        self.x              = x
        self.y              = y
        self.y_error        = y_error
        self.y_error_scaled = 0
        self.red_chi_squ    = 1
        self.r_squ          = 1
        self.dof            = 0
    
    def get_unweighted_model(self, start):
        p_opt, p_cov = curve_fit(self.function, xdata=self.x, ydata=self.y, 
                             p0=start, absolute_sigma=False, maxfev=10000)
    
        p_err = np.sqrt(np.diag(p_cov))

        fit_a = ufloat(p_opt[0], p_err[0])
        fit_b = ufloat(p_opt[1], p_err[1])
        fit_c = ufloat(p_opt[2], p_err[2])
    
        ## Calculation of R^2
        residuals = self.y - self.function(self.x, *p_opt) 
        ss_res = np.sum(residuals**2) 
        ss_tot = np.sum((self.y - np.mean(self.y))**2) 
        
        self.r_squ = 1 - (ss_res / ss_tot)
        self.dof   = len(self.x) - 3
        
        errors_scaled = np.sum(((self.y - self.function(self.x, *p_opt)))**2)**0.5
        self.y_error_scaled = unumpy.std_devs(unumpy.uarray(unumpy.nominal_values(self.y), errors_scaled))
        
        print("Unweighted fit with three parameters:")
        print("1st parameter =", fit_a, " unit")
        print("2nd parameter =", fit_b, " unit")
        print("3rd parameter =", fit_c, " unit")
        print("R^2 (2 parameter model) =", round(self.r_squ,4), ", dof =", self.dof)
        print(" ")
    
        return fit_a, fit_b, fit_c

    def unweighted(self, start):
        fit_a, fit_b, fit_c = self.get_unweighted_model(start)
    
        x_fit = np.linspace(self.x[0], self.x[len(self.x)-1], 1000)
        y_fit = self.function(x_fit, fit_a.n, fit_b.n, fit_c.n)
    
        return np.array([x_fit, y_fit]), fit_a, fit_b, fit_c


    def get_weighted_model(self, start):
        p_opt, p_cov = curve_fit(self.function, xdata=self.x, ydata=self.y, sigma=self.y_error, 
                                 p0=start, absolute_sigma=False, maxfev=10000)
    
        p_err = np.sqrt(np.diag(p_cov))

        fit_a = ufloat(p_opt[0], p_err[0])
        fit_b = ufloat(p_opt[1], p_err[1])
        fit_c = ufloat(p_opt[2], p_err[2])
    
        ## Calculation of R^2
        residuals = self.y - self.function(self.x, *p_opt) 
        ss_res = np.sum(residuals**2) 
        ss_tot = np.sum((self.y - np.mean(self.y))**2) 

        self.r_squ = 1 - (ss_res / ss_tot)
        self.dof   = len(self.x) - 3
        
        self.red_chi_squ = np.sum(((self.y - self.function(self.x, *p_opt)) / self.y_error)**2) / self.dof
        
        errors_scaled = np.sum(((self.y - self.function(self.x, *p_opt)))**2)**0.5
        self.y_error_scaled = unumpy.std_devs(unumpy.uarray(unumpy.nominal_values(self.y), errors_scaled))
        
        print("Weighted fit with three parameters:")
        print("1st parameter =", fit_a, " unit")
        print("2nd parameter =", fit_b, " unit")
        print("3rd parameter =", fit_c, " unit")
        print("R^2 (2 parameter model) =", round(self.r_squ,4), ", dof = ", self.dof)
        print("red Chi^2 (2 parameter model) =", round(self.red_chi_squ,4), ", dof = ", self.dof)
        print(" ")
    
        return fit_a, fit_b, fit_c

    def weighted_fit_params_three(self, start):
        fit_a, fit_b, fit_c = self.get_weighted_model(start)
    
        x_fit = np.linspace(self.x[0], self.x[len(self.x)-1], 1000)
        y_fit = self.function(x_fit, fit_a.n, fit_b.n, fit_c.n)
    
        return np.array([x_fit, y_fit]), fit_a, fit_b, fit_c

# ┌──────────────────────┐
# │ Export data to Latex │
# └──────────────────────┘
class export_latex:
    def __init__(self, pd_data_frame, file_name):
        self.pd_data_frame = pd_data_frame
        self.file_name     = file_name
        
    def save_dataframe(self):
        text_file = open(path+"/tables/"+self.file_name+".txt", "w")
        nothing_important = text_file.write(self.pd_data_frame.to_latex(index=False))
        text_file.close()
