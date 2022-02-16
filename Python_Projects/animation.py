#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 23:12:18 2020

@author: Florian Kluibenschedl
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib
#matplotlib.use("Agg")

import glob
import os, sys

path = os.path.dirname(__file__)
sys.path.append(path)

import calculations as calc # for more involved calculations

# ┌───────────┐
# │ Read data │
# └───────────┘
input_object = calc.handle_input()
folder_path  = input_object.get_params()

xData     = []
yData     = []
thetaData = []
pause = False

t = np.linspace(0, int(input_object.dt*input_object.Nsim), int(input_object.Nsim/input_object.Nsave)+1)
for i in t:
    x     = np.array(pd.read_csv(folder_path + "/configuration_t_{0}00000".format(i), usecols=[0], delimiter=" "))[:, -1]
    y     = np.array(pd.read_csv(folder_path + "/configuration_t_{0}00000".format(i), usecols=[1], delimiter=" "))[:, -1]
    theta = np.array(pd.read_csv(folder_path + "/configuration_t_{0}00000".format(i), usecols=[2], delimiter=" "))[:, -1]

    xData.append(x)
    yData.append(y)
    thetaData.append(theta)

intervall_rate = input_object.get_intervall_rate()

# ┌─────────┐
# │ Animate │
# └─────────┘
fig = plt.figure()
ax = fig.add_subplot(111)

def onClick(event):
    global pause
    pause ^= True

def update(frame):
    if not pause:
        ax.cla()
        x_data     = xData[frame]
        y_data     = yData[frame]
        theta_data = thetaData[frame]

        ax.quiver(x_data, y_data, np.cos(theta_data), np.sin(theta_data))
        ax.set_title("t = {0}".format(t[frame]))

        ax.set_xticks([])
        ax.set_yticks([])

fig.canvas.mpl_connect('button_press_event', onClick)
animator = ani.FuncAnimation(fig, update, frames=range(len(xData)), interval=intervall_rate)

file_name = folder_path + r"/animation.mp4"
writervideo = ani.FFMpegWriter(fps=60)
animator.save(file_name, writer=writervideo)

plt.show()
plt.close()
