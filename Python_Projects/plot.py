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
import matplotlib.pyplot as plt

# ┌─────────────┐
# │ Without Fit │
# └─────────────┘
class without_fit_three_data_lines:
    def __init__(self, x, y1, y2, y3, y_error, x_label, y_label, data_label):
        self.x          = x
        self.y1         = y1
        self.y2         = y2
        self.y3         = y3
        self.y_error    = y_error
        self.x_label    = x_label
        self.y_label    = y_label
        self.data_label = data_label

    def scatter(self, image_name, set_grid, set_legend):
        plt.plot(self.x, self.y1, label=r"$E_{pot}$")
        plt.plot(self.x, self.y2, label=r"$E_{kin}$")
        plt.plot(self.x, self.y3, label=r"$E_{sum}$")

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(path+"/cs07_files/"+image_name+".png", dpi=800)
        plt.show()

    def scatter_error(self, image_name, set_grid, set_legend):
        plt.errorbar(self.x, self.y, yerr=self.y_error, capsize=3, marker='.', linestyle = 'None', label=self.data_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(path+"/images/"+image_name+".png", dpi=800)
        plt.show()

class phase_analysis:
    def __init__(self, x, y, params, x_label, y_label, plot_label):
        self.x          = x
        self.y          = y
        self.params     = params
        self.x_label    = x_label
        self.y_label    = y_label
        self.plot_label = plot_label

    def plot_phases(self, param_label, image_name, set_grid, set_legend):
        plt.title(self.plot_label)
        for i in range(len(self.params)):
            plt.plot(self.x, self.y[i], label=param_label + r"$ = $" + str(self.params[i]))

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        plt.ylim(ymin=0, ymax=1.1)
        plt.axhline(y=1, color="grey", linestyle="dashed", linewidth=1)

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend(loc="upper right")

        plt.savefig(image_name+".png", dpi=800)
        plt.show()

class without_fit_one_data_line:
    def __init__(self, x, y, x_label, y_label, data_label):
        self.x          = x
        self.y          = y
        self.x_label    = x_label
        self.y_label    = y_label
        self.data_label = data_label

    def scatter(self, image_name, set_grid, set_legend):
        plt.scatter(self.x, self.y, marker="x", label=self.data_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(image_name+".png", dpi=800)
        plt.show()

    def plot_phase(self, image_name, set_grid, set_legend):
        fig= plt.figure()
        ax = plt.subplot()

        plt.plot(self.x, self.y, label=self.data_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        plt.ylim(ymin=0, ymax=1.1)
        plt.axhline(y=1, color="grey", linestyle="dashed", linewidth=1)

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(image_name+".png", dpi=800)
        plt.show()

    def no_scatter(self, image_name, set_grid, set_legend):
        plt.plot(self.x, self.y, label=self.data_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(image_name+".png", dpi=800)
        plt.show()

    def scatter_error(self, y_error, image_name, set_grid, set_legend):
        plt.errorbar(self.x, self.y, yerr=y_error, capsize=3, marker='.', linestyle = 'None', label=self.data_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(image_name+".png", dpi=800)
        plt.show()

# ┌──────────┐
# │ With Fit │
# └──────────┘
class with_fit:
    def __init__(self, x, y, y_error, x_fit, y_fit, x_label, y_label, data_label, fit_label):
        self.x          = x
        self.y          = y
        self.y_error    = y_error
        self.x_fit      = x_fit
        self.y_fit      = y_fit
        self.x_label    = x_label
        self.y_label    = y_label
        self.data_label = data_label
        self.fit_label  = fit_label

    def scatter_unweighted(self, image_name, set_grid, set_legend):
        plt.scatter(self.x, self.y, color="black", label=self.data_label)
        plt.plot(self.x_fit, self.y_fit, color="red", label=self.fit_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(path+"/images/"+image_name+".png", dpi=800)
        plt.show()

    def scatter_weighted(self, image_name, set_grid, set_legend):
        plt.errorbar(self.x, self.y, yerr=self.y_error, capsize=3, marker='.', linestyle = 'None', label=self.data_label)
        plt.plot(self.x_fit, self.y_fit, color="red", label=self.fit_label)

        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        if set_grid:
            plt.grid()
        if set_legend:
            plt.legend()

        plt.savefig(path+"/images/"+image_name+".png", dpi=800)
        plt.show()

# ┌───────────────┐
# │ Configuration │
# └───────────────┘
class configuration:
    def __init__(self, x, y, theta, x_label, y_label, data_label):
        self.x          = x
        self.y          = y
        self.theta      = theta
        self.x_label    = x_label
        self.y_label    = y_label
        self.data_label = data_label

    def plot(self, image_name, set_grid, set_legend):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.quiver(self.x, self.y, np.cos(self.theta), np.sin(self.theta))
        #ax.set_title("t = {0}".format(t[frame]))

        plt.savefig(image_name+".png", dpi=800)
        plt.show()
