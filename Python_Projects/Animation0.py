import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as ani
import os
import glob
import matplotlib
#matplotlib.use("Agg")

xData     = []
yData     = []
thetaData = []
pause = False

# assign directory
directory = '../simulation_results/N_1000_L_1.000000_v_0.030000_R_1.000000_D_0.100000/Neq_2000_Nsim_1000_dt_1.000000/'

t = np.linspace(1000, 100000, 100)
t = np.linspace(10, 1000, 100)

for i in t:
    x     = np.array(pd.read_csv(directory + "configuration_t_{0}00000".format(i), usecols=[0], delimiter=" "))[:, -1]
    y     = np.array(pd.read_csv(directory + "configuration_t_{0}00000".format(i), usecols=[1], delimiter=" "))[:, -1]
    theta = np.array(pd.read_csv(directory + "configuration_t_{0}00000".format(i), usecols=[2], delimiter=" "))[:, -1]

    xData.append(x)
    yData.append(y)
    thetaData.append(theta)

o = np.array([[0, 0], [0, 0]])

X = np.array((0, 0))
Y = np.array((0, 0))
U = np.array([1, 2])
V = np.array([-2, 4])

# fig, ax = plt.subplots(figsize=(16, 9))
# q = ax.quiver(test[0], test[1], test[2], test[3], units='xy', scale=1, width=0.05)

# print(test)

# ax.set_aspect('equal')

# plt.xlim(0, 20)
# plt.ylim(0, 20)
# plt.show()

# fig, ax = plt.subplots(figsize=(16,9))
# plt.plot(xData[0], yData[0], 'k.')
# plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)


def onClick(event):
    global pause
    pause ^= True


def update(frame):
    if not pause:
        ax.cla()
        # ax.plot(xData[frame], yData[frame], "k.")
        x_data     = xData[frame]
        y_data     = yData[frame]
        theta_data = thetaData[frame]

        ax.quiver(x_data, y_data, np.cos(theta_data), np.sin(theta_data))
        ax.set_title("t = {0}".format(t[frame]))

fig.canvas.mpl_connect('button_press_event', onClick)
animator = ani.FuncAnimation(fig, update, frames=range(len(xData)), interval=10)

f = r"animation.mp4"
writervideo = ani.FFMpegWriter(fps=60)
animator.save(f, writer=writervideo)

plt.show()
plt.close()
