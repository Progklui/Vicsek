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
        self.folder_pa = []
        self.N_1       = 1
        self.N_2       = 1
        self.N_inc     = 1
        self.L         = 1
        self.dt        = 1
        self.R_1       = 1
        self.R_2       = 1
        self.R_inc     = 1
        self.v_1       = 1
        self.v_2       = 1
        self.v_inc     = 1
        self.D_rot_1   = 1
        self.D_rot_2   = 1
        self.D_rot_inc = 1
        self.Neq       = 1
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

        self.N_1       = params[0]
        self.N_2       = params[1]
        self.N_inc     = params[2]
        self.L         = params[3]
        self.dt        = params[4]
        self.R_1       = params[5]
        self.R_2       = params[6]
        self.R_inc     = params[7]
        self.v_1       = params[8]
        self.v_2       = params[9]
        self.v_inc     = params[10]
        self.D_rot_1   = params[11]
        self.D_rot_2   = params[12]
        self.D_rot_inc = params[13]
        self.Nsim      = params[14]
        self.Nsave     = params[15]
        self.Neq       = params[16]

        # sim_params = "N_" + str(round(self.N)) + "_L_" + '{:.6f}'.format(self.L) + "_v_" + '{:.6f}'.format(self.v) + "_R_" + '{:.6f}'.format(self.R) + "_D_" + '{:.6f}'.format(self.D_rot)
        # sim_counts = "Neq_" + str(2000) + "_Nsim_" + str(round(self.Nsim)) + "_dt_" + '{:.6f}'.format(self.dt)

        # self.folder_pa = "../simulation_results/" + sim_params + "/" + sim_counts

        print(" ")
        print("Verify input:")
        print("Parameter file:", "../Parameter/"+self.file_path)
        print("Folder path  : ", self.folder_pa)
        print(" ")
        print("N     =        ", self.N_1, "-", self.N_2, ", step size: ", self.N_inc)
        print("L     =        ", self.L)
        print("dt    =        ", self.dt)
        print("R     =        ", self.R_1, "-", self.R_2, ", step size: ", self.R_inc)
        print("v     =        ", self.v_1, "-", self.v_2, ", step size: ", self.v_inc)
        print("D_rot =        ", self.D_rot_1, "-", self.D_rot_2, ", step size: ", self.D_rot_inc)
        print("Neq   =        ", self.Neq)
        print("Nsim  =        ", self.Nsim)
        print("Nsave =        ", self.Nsave)
        print(" ")
        accept = input("Analyze all (y/n/break)? ")
        if accept == "y" or accept == "yes":
            print(" ")
            folder_path = []

            N_array     = np.arange(self.N_1, self.N_2, self.N_inc)
            v_array     = np.arange(self.v_1, self.v_2, self.v_inc)
            R_array     = np.arange(self.R_1, self.R_2, self.R_inc)
            D_rot_array = np.arange(self.D_rot_1, self.D_rot_2, self.D_rot_inc)

            for N in N_array:
                for v in v_array:
                    for R in R_array:
                        for D_rot in D_rot_array:
                            folder_path.append(self.get_folder_path(N, self.L, v, R, D_rot, self.Neq, self.Nsim, self.dt))
            return folder_path
        elif accept == "n" or accept == "no":
            print(" ")
            N     = int(input("N     = "))
            L     = float(input("L     = "))
            dt    = float(input("dt    = "))
            R     = float(input("R     = "))
            v     = float(input("v     = "))
            D_rot = float(input("D_rot = "))
            Neq   = int(input("Neq   = "))
            Nsim  = int(input("Nsim  = "))
            print(" ")

            folder_path = []
            folder_path.append(self.get_folder_path(N, L, v, R, D_rot, Neq, Nsim, dt))
            return folder_path
        else:
            exit()

    def get_folder_path(self, N, L, v, R, D_rot, Neq, Nsim, dt):
        sim_params = "N_" + str(round(N)) + "_L_" + '{:.6f}'.format(L) + "_v_" + '{:.6f}'.format(v) + "_R_" + '{:.6f}'.format(R) + "_D_" + '{:.6f}'.format(D_rot)
        sim_counts = "Neq_" + str(round(Neq)) + "_Nsim_" + str(round(Nsim)) + "_dt_" + '{:.6f}'.format(dt)

        return "../simulation_results/" + sim_params + "/" + sim_counts

    def get_configuration(self):
        accept = input("Print (other) configuration (y/n)? ")
        if accept == "y" or accept == "yes":
            time = input("Specify time: ")
            file_name = "/configuration_t_" + time
            print(" ")
            return file_name, time, True
        elif accept == "n" or accept == "no":
            print(" ")
            return "empty", 12, False
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
