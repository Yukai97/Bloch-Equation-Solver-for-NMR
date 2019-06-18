import numpy as np
from numpy import pi
from matplotlib import pyplot as plt
from Bloch_Equation_Solver.Bloch_equations_solver import BlochEquationsSolver


def gaussian(t, t0, sigma, A):
    y = A*np.exp(-(t - t0)**2/(2*sigma**2))
    return y


# parameters
gamma_He = 2*pi*32.43409966   # MHz/T
gamma_Ne = gamma_He/9.649
gamma_Xe = 2*pi*11.77717
B0 = 2*pi*15.0148/gamma_He
omega_He = gamma_He*B0
omega_Ne = gamma_Ne*B0
omega_Xe = gamma_Xe*B0

# Adjustable parameters
Ne_BxA = 0.02
Xe_BxA = 0.015
He_BxA = 0.005
gaussianA = 1
t0_Ne = 6.3
t0_Xe = 0
t0_He = 0
sigma_Ne = 0.91
sigma_Xe = 2.26
sigma_He = 2.468

# Initial Magnetization
Ne_Mx0 = 0
Ne_My0 = 0
Ne_Mz0 = 1
Xe_Mx0 = 0
Xe_My0 = 0
Xe_Mz0 = 1
He_Mx0 = 0
He_My0 = 0
He_Mz0 = 1

# Build time array. When solving tipping pulse we can use very small time_step to minimize numerical errors.
t_total = 10
time_step = 0.0002
time_array = np.arange(0, t_total + time_step, time_step)

# Magnetic field
index_ne = int(t0_Ne/time_step)
bx_ne1 = Ne_BxA * np.sin(omega_Ne * time_array[0:index_ne])
bx_ne2 = Ne_BxA * np.sin(omega_Ne * time_array[index_ne:]) * gaussian(time_array[index_ne:], t0_Ne, sigma_Ne, gaussianA)
bx_ne = np.concatenate([bx_ne1, bx_ne2], axis=None)
bx_xe = Xe_BxA * np.sin(omega_Xe * time_array) * gaussian(time_array, t0_Xe, sigma_Xe, gaussianA)
bx_he = He_BxA * np.sin(omega_He * time_array) * gaussian(time_array, t0_He, sigma_He, gaussianA)
Bx = bx_ne + bx_he + bx_xe

By = np.zeros(len(time_array))
Bz = B0 * np.ones(len(time_array))

ne_bes = BlochEquationsSolver(gamma_Ne, time_array, time_step, [Ne_Mx0, Ne_My0, Ne_Mz0], [Bx, By, Bz])
he_bes = BlochEquationsSolver(gamma_He, time_array, time_step, [He_Mx0, He_My0, He_Mz0], [Bx, By, Bz])
xe_bes = BlochEquationsSolver(gamma_Xe, time_array, time_step, [Xe_Mx0, Xe_My0, Xe_Mz0], [Bx, By, Bz])
plt.figure()
plt.plot(time_array, ne_bes.M[2], time_array, xe_bes.M[2], time_array, he_bes.M[2])
plt.show()
