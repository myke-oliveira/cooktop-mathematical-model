#!/usr/bin/env python
# coding: utf-8

# Espira em Coordenadas Polares
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from mpl_toolkits.mplot3d import Axes3D
from time import sleep

# Parameters
R = 5e-2
N = 10
e = 5e-3
d = 1e-2
numsteps = 300

# Variables
t = np.linspace(0, 1)
omega = 2 * np.pi
i = np.sin(omega*t)
#i = 1

# r_L, theta: coordenadas polares no condutor
theta = np.linspace(0, 2*N*np.pi, numsteps)
r_L = R - e/(2*np.pi)*theta

fig = plt.figure(1)
plt.polar(theta, r_L)
plt.grid()
fig.show(False)
sleep(.5)

# Espira em Coordenadas Cartesianas

# x_L, y_L, z_L: coordenadas retângulares no condutor
x_L = r_L * np.cos(theta)
y_L = r_L * np.sin(theta)
z_L = np.zeros(len(x_L))

fig = plt.figure(2)
plt.plot(x_L, y_L)
fig.show(False)
sleep(.5)

# # Elementos do condutor orientado

#dx_L, dy_L dz_L: coordenadas dos elementos do contudor
dX_L = x_L[1:] - x_L[:-1]
dY_L = y_L[1:] - y_L[:-1]
dZ_L = np.zeros(len(dX_L))

# X_L, Y_L: coordenadas dos pontos médios dos elementos do condutor
X_L = (x_L[1:] + x_L[:-1])/2
Y_L = (y_L[1:] + y_L[:-1])/2
Z_L = np.zeros(len(X_L))

fig, ax = plt.subplots()
q = ax.quiver(X_L, Y_L, dX_L, dY_L)
fig.show(False)


for tnow, inow in zip(t, i):
    x_Blist = list()
    y_Blist = list()
    z_Blist = list()
    x_Plist = list()
    y_Plist = list()
    z_Plist = list()
    for x_P, y_P, z_P in product(np.arange(-.06, .06, .01), np.arange(-.06, .06, .01), np.arange(0, .06, .01)):
        x_Plist.append(x_P)
        y_Plist.append(y_P)
        z_Plist.append(z_P)
        X_R = (x_P- x_0 for x_0 in x_L)
        Y_R = (y_P - y_0 for y_0 in y_L)
        Z_R = (y_P - z_0 for z_0 in z_L)
        # print(*X_R, *Y_R, *Z_R)
        B = 0
        for x_r, y_r, z_r, dx_L, dy_L, dz_L in zip(X_R, Y_R, Z_R, dX_L, dY_L, dZ_L):
            r_squared = x_r**2 + y_r**2 + z_r**2
            x_r_hat = x_r / r_squared**5
            y_r_hat = y_r / r_squared**5
            z_r_hat = z_r / r_squared**5
            dl = np.array((dx_L, dy_L, dz_L))
            r_hat = np.array((x_r_hat, y_r_hat, z_r_hat))
            B += 1e-7*inow*np.cross(dl, r_hat)/r_squared
        x_Blist.append(B[0])
        y_Blist.append(B[1])
        z_Blist.append(B[2])
    fig = plt.figure()
    aux = fig.gca(projection='3d')
    aux.quiver(x_Plist, y_Plist, z_Plist, x_Blist, y_Blist, z_Blist, length=.01, normalize=True)
    aux.quiver(X_L, Y_L, Z_L, dX_L, dY_L, dZ_L, length=.01, normalize=True, color='yellow')
    plt.savefig(f'figure {tnow}.png')
input()

