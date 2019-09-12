from sympy import *
from scipy.integrate import quad
import numpy as np

init_printing()

R, N, e, d, i = symbols('R N e d i(t)')
t, x, y, z, = symbols('t x y z')
r, theta = symbols('r theta')
E = Symbol('E(x, y, z, t)')
B = Symbol('B(x, y, z, t)')
J = Symbol('J(x, y, z, t)')
P = Symbol('P(x, y, z, t)')
epsilon, mu = symbols('epsilon mu')

r_L = R - e/(2*pi) * theta
print('-'*80)
print('r_l = ')
print(r_L)

x_L = r_L * cos(theta)
y_L = r_L * sin(theta)
print('-'*80)
print('x_L = ')
print(x_L)
print('-'*80)
print('y_L = ')
print(y_L)

r1 = x - x_L
r2 = y - y_L
r3 = z
print('-'*80)
print('r1 = ')
print(r1)
print('r2 = ')
print(r2)
print('r3 = ')
print(r3)

r_squared = r1**2 + r2**2 + r3**2
print('r**2 = ')
print(r_squared)

r_mod = sqrt(r_squared)
print('|r| = ')
print(r_mod)

r1_hat = r1 / r_mod
r2_hat = r2 / r_mod
r3_hat = r3 / r_mod
print('r1_hat = ')
print(r1_hat)
print('r2_hat = ')
print(r2_hat)
print('r3_hat = ')
print(r3_hat)

dl1_dtheta = diff(x_L, theta)
dl2_dtheta = diff(y_L, theta)
dl3_dthetha = 0
print('dl1/dtheta = ')
print(dl1_dtheta)
print('dl2/dtheta = ')
print(dl2_dtheta)
print('dl3/dtheta = ')
print(dl3_dthetha)

c = mu / (2*pi)
print('c = ')
print(c)

dl_times_r_over_r_squared_1 = dl2_dtheta * r3_hat - r2_hat * dl3_dthetha
dl_times_r_over_r_squared_2 = -dl1_dtheta * r3_hat - r1_hat * dl3_dthetha
dl_times_r_over_r_squared_3 = dl1_dtheta * r2_hat - r1_hat * dl2_dtheta

dB1_dtheta = c*dl_times_r_over_r_squared_1
dB2_dtheta = c*dl_times_r_over_r_squared_2
dB3_dtheta = c*dl_times_r_over_r_squared_3

print('dB1/dtheta = ')
print(dB1_dtheta)
print('dB2/dtheta = ')
print(dB2_dtheta)
print('dB3/dtheta = ')
print(dB3_dtheta)

# int_dl_times_r_over_r_squared_1 = integrate(dl_times_r_over_r_squared_1, theta)
# int_dl_times_r_over_r_squared_2 = integrate(dl_times_r_over_r_squared_2, theta)
# int_dl_times_r_over_r_squared_3 = integrate(dl_times_r_over_r_squared_3, theta)

print('Numeric Integration at (0, 0, 0)')
#R = float(input('R= '))
#N = float(input('N= '))
#e = float(input('e = '))
#d = float(input('d = '))
#i = float(input('i = '))

newR, newN, newe, newd, newi = 5, 10, .5, .1, 1
newmu, newepsilon = 4e-7*pi, 8.85e-12

dB1_dtheta = dB1_dtheta.subs(R, newR)
dB1_dtheta = dB1_dtheta.subs(N, newN)
dB1_dtheta = dB1_dtheta.subs(e, newe)
dB1_dtheta = dB1_dtheta.subs(z, newd)
dB1_dtheta = dB1_dtheta.subs(i, newi)
dB1_dtheta = dB1_dtheta.subs(epsilon, newepsilon)
dB1_dtheta = dB1_dtheta.subs(mu, newmu)
dB1_dtheta = dB1_dtheta.subs(x, newx)
dB1_dtheta = dB1_dtheta.subs(y, newy)
dB1_dtheta = dB1_dtheta.subs(z, newz)
dB2_dtheta = dB2_dtheta.subs(R, newR)
dB2_dtheta = dB2_dtheta.subs(N, newN)
dB2_dtheta = dB2_dtheta.subs(e, newe)
dB2_dtheta = dB2_dtheta.subs(z, newd)
dB2_dtheta = dB2_dtheta.subs(i, newi)
dB2_dtheta = dB2_dtheta.subs(epsilon, newepsilon)
dB2_dtheta = dB2_dtheta.subs(mu, newmu)
dB2_dtheta = dB2_dtheta.subs(x, newx)
dB2_dtheta = dB2_dtheta.subs(y, newy)
dB2_dtheta = dB2_dtheta.subs(z, newz)
dB3_dtheta = dB3_dtheta.subs(R, newR)
dB3_dtheta = dB3_dtheta.subs(N, newN)
dB3_dtheta = dB3_dtheta.subs(e, newe)
dB3_dtheta = dB3_dtheta.subs(z, newd)
dB3_dtheta = dB3_dtheta.subs(i, newi)
dB3_dtheta = dB3_dtheta.subs(epsilon, newepsilon)
dB3_dtheta = dB3_dtheta.subs(mu, newmu)
dB3_dtheta = dB3_dtheta.subs(x, newx)
dB3_dtheta = dB3_dtheta.subs(y, newy)
dB3_dtheta = dB3_dtheta.subs(z, newz)

print('dB1/dtheta = ')
print(dB1_dtheta)
print('dB2/dtheta = ')
print(dB2_dtheta)
print('dB3/dtheta = ')
print(dB3_dtheta)

B1 = quad(lambda x: dB1_dtheta.subs(theta, x).evalf(), 0, 2*newN*pi.evalf())
print('B1 = ')
print(B1)
B2 = quad(lambda x: dB2_dtheta.subs(theta, x).evalf(), 0, 2*newN*pi.evalf())
print('B2 = ')
print(B2)
B3 = quad(lambda x: dB3_dtheta.subs(theta, x).evalf(), 0, 2*newN*pi.evalf())
print('B3 = ')
print(B3)

