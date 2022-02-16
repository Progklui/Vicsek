#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 11:58:35 2020

@author: Florian Kluibenschedl
"""
import numpy as np
import pandas as pd

import os, sys
from sys import exit

path = os.path.dirname(__file__)
sys.path.append(path)

class handle_input:
    def __init__(self):
        self.file_path = "input.csv"
        self.folder_pa = "N_"
        self.N         = 1
        self.L         = 1
        self.dt        = 1
        self.R         = 1
        self.v         = 1
        self.D_rot     = 1
        self.Nsim      = 1
        self.Nsave     = 1

    def help_function(self):
        print(" ")
        print("Documentation of parameters:")
        print("Argument structure: [PATH]")
        print(" ")
        print("[PATH] = ../Parameter/input.csv (=DEFAULT)")
        print(" ")
        quit()

    def get_params(self):
        try:
            argument = sys.argv[1]
            if argument == "-h":
                self.help_function()
            else:
                self.file_path = argument
        except:
            pass

        params     = np.array(pd.read_csv("../Parameter/"+self.file_path, usecols=[1], delimiter=";"))[:, -1] # np.loadtxt(self.file_path, delimiter=';', skiprows=1).T

        self.N     = params[0]
        self.L     = params[1]
        self.dt    = params[2]
        self.R     = params[3]
        self.v     = params[4]
        self.D_rot = params[5]
        self.Nsim  = params[6]
        self.Nsave = params[7]

        sim_params = "N_" + str(round(self.N)) + "_L_" + '{:.6f}'.format(self.L) + "_v_" + '{:.6f}'.format(self.v) + "_R_" + '{:.6f}'.format(self.R) + "_D_" + '{:.6f}'.format(self.D_rot)
        sim_counts = "Neq_" + str(2000) + "_Nsim_" + str(round(self.Nsim)) + "_dt_" + '{:.6f}'.format(self.dt)

        self.folder_pa = "../simulation_results/" + sim_params + "/" + sim_counts

        print(" ")
        print("Verify input:")
        print("Paramter file: ", "../Parameter/"+self.file_path)
        print("Folder path  : ", self.folder_pa)
        print(" ")
        print("N     =        ", self.N)
        print("L     =        ", self.L)
        print("dt    =        ", self.dt)
        print("R     =        ", self.R)
        print("v     =        ", self.v)
        print("D_rot =        ", self.D_rot)
        print("Nsim  =        ", self.Nsim)
        print("Nsave =        ", self.Nsave)
        print(" ")
        accept = input("Accept input (y/n)? ")
        if accept == "y" or accept == "yes":
            print(" ")
            return self.folder_pa
        else:
            exit()

    def get_configuration(self):
        accept = input("Print (other) configuration (y/n)? ")
        if accept == "y" or accept == "yes":
            time = input("Specify time: ")
            file_name = "/configuration_t_" + time
            print(" ")
            return file_name, time, True
        else:
            exit()

    def get_trajectory(self):
        accept = input("Print trajectory of order parameter (y/n)? ")
        if accept == "y" or accept == "yes":
            print(" ")
            return True
        elif accept == "n" or accept == "no":
            print(" ")
            return False
        else:
            exit()

    def get_intervall_rate(self):
        intervall_rate = input("Intervall = ")
        return int(intervall_rate)
