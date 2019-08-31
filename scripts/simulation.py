#!/usr/bin/env python
# coding: utf-8

# Espira em Coordenadas Polares
import numpy as np
import matplotlib.pyplot as plt

# Parameters
R = 5e-2
N = 10
e = 5e-3
d = 1e-2
numsteps = 300

# Variables
i = 1.

# r_L, theta: coordenadas polares no condutor
theta = np.linspace(0, 2*N*np.pi, numsteps)
r_L = R - e/(2*np.pi)*theta

plt.polar(theta, r_L)
plt.grid()
plt.show(True)


# Espira em Coordenadas Cartesianas

# x_L, y_L, z_L: coordenadas retângulares no condutor
x_L = r_L * np.cos(theta)
y_L = r_L * np.sin(theta)
z_L = np.zeros(len(x_L))

plt.plot(x_L, y_L)
plt.show(True)

# # Elementos do condutor orientado

#dx_L, dy_L dz_L: coordenadas dos elementos do contudor
dX_L = x_L[1:] - x_L[:-1]
dY_L = y_L[1:] - y_L[:-1]
dZ_L = np.zeros(len(dX_L))

# X_L, Y_L: coordenadas dos pontos médios dos elementos do condutor
X_L = (x_L[1:] + x_L[:-1])/2
Y_L = (y_L[1:] + y_L[:-1])/2

fig, ax = plt.subplots()
q = ax.quiver(X_L, Y_L, dX_L, dY_L)
fig.show(True)

x_P, y_P, z_P = 3e-2, 3e-2, 3e-2
X_R = (x_P- x_0 for x_0 in x_L)
Y_R = (y_P - y_0 for y_0 in y_L)
Z_R = (y_P - z_0 for z_0 in z_L)
# print(*X_R, *Y_R, *Z_R)

B = 0
#for x_r, y_r, z_r, dx_L, dy_L, dz_L in zip(X_R, Y_R, Z_R, dX_L, dY_L, dZ_L):
#    r_squared = x_r**2 + y_r**2 + z_r**2
#    x_r_hat = x_r / r_squared**5
#    y_r_hat = y_r / r_squared**5
#    z_r_hat = z_r / r_squared**5
#    dl = np.array((dx_L, dy_L, dz_L))
#    r_hat = np.array((x_r_hat, y_r_hat, z_r_hat))
    
#    

