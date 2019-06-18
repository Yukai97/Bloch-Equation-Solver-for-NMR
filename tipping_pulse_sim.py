from Bloch_Equation_Solver.Bloch_equations_solver import BlochEquationsSolver
import numpy as np
import math


def gaussian(t, t0, sigma, A):
    y = A*math.exp(-(t - t0)**2/(2*sigma**2))
    return y


class TippingPulseSimulator:

    def __init__(self, he_angle, ne_angle, M0_he, M0_ne, gauss_A=1):
        self.he_angle = he_angle
        self.ne_angle = ne_angle
        self.xe_angle = 90
        self.gauss_A = gauss_A
        self.gamma = {'H': 2*math.pi*32.43409966, 'X': 2*math.pi*11.77717}
        self.gamma['N'] = self.gamma['H']/9.649
        self.B0 = 2*math.pi*15.0148/self.gamma['H']
        self.omega = {'H': self.gamma['H'] * self.B0, 'N': self.gamma['N'] * self.B0, 'X': self.gamma['X'] * self.B0}
        self.time_step = 0.001

        self.angle_pattern = self.set_angle_pattern()
        self.look_up_table = self.build_look_up_table()
        self.pulse_paras = self.look_up_table[self.angle_pattern]
        time_array, pulse_B = self.generate_tipping_pulse()
        self.time_array = time_array
        self.pulse_B = pulse_B
        By, Bz = self.set_magnetic_field()

        self.he_bes = BlochEquationsSolver(self.gamma['H'], time_array, self.time_step, M0_he, [pulse_B, By, Bz])
        self.ne_bes = BlochEquationsSolver(self.gamma['N'], time_array, self.time_step, M0_ne, [pulse_B, By, Bz])

    def set_angle_pattern(self):
        angle_pattern = 'H' + str(self.he_angle) + 'N' + str(self.ne_angle) + 'X' + str(self.xe_angle)
        return angle_pattern

    def build_look_up_table(self):
        look_up_table = {}
        look_up_table['H20N90X90'] = {'t0_ne': 6.3, 'sigma_he': 0.55, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H40N90X90'] = {'t0_ne': 6.3, 'sigma_he': 1.093, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H60N90X90'] = {'t0_ne': 6.3, 'sigma_he': 1.642, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H80N90X90'] = {'t0_ne': 6.3, 'sigma_he': 2.192, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H90N90X90'] = {'t0_ne': 6.3, 'sigma_he': 2.055, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                      'B_he': 0.006, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H100N90X90'] = {'t0_ne': 6.3, 'sigma_he': 2.286, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                       'B_he': 0.006, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H120N90X90'] = {'t0_ne': 6.3, 'sigma_he': 2.75, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                       'B_he': 0.006, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H140N90X90'] = {'t0_ne': 6.3, 'sigma_he': 3.23, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                       'B_he': 0.006, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H160N90X90'] = {'t0_ne': 6.3, 'sigma_he': 3.195, 'sigma_ne': 0.91, 'sigma_xe': 2.26,
                                       'B_he': 0.007, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H90N20X90'] = {'t0_ne': 0, 'sigma_he': 2.46, 'sigma_ne': 1.32, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H90N40X90'] = {'t0_ne': 2, 'sigma_he': 2.465, 'sigma_ne': 1.04, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H90N60X90'] = {'t0_ne': 3.5, 'sigma_he': 2.467, 'sigma_ne': 1.16, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H90N80X90'] = {'t0_ne': 5, 'sigma_he': 2.468, 'sigma_ne': 1.285, 'sigma_xe': 2.26,
                                      'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 10}
        look_up_table['H90N100X90'] = {'t0_ne': 7, 'sigma_he': 2.468, 'sigma_ne': 1.01, 'sigma_xe': 2.26,
                                       'B_he': 0.005, 'B_ne': 0.02, 'B_xe': 0.015, 't_total': 12}
        look_up_table['H90N120X90'] = {'t0_ne': 8, 'sigma_he': 2.47, 'sigma_ne': 1.155, 'sigma_xe': 2.26,
                                       'B_he': 0.005, 'B_ne': 0.021, 'B_xe': 0.015, 't_total': 12}
        look_up_table['H90N140X90'] = {'t0_ne': 7, 'sigma_he': 2.483, 'sigma_ne': 1.52, 'sigma_xe': 2.263,
                                       'B_he': 0.005, 'B_ne': 0.026, 'B_xe': 0.015, 't_total': 12}
        look_up_table['H90N160X90'] = {'t0_ne': 6.6, 'sigma_he': 2.512, 'sigma_ne': 1.33, 'sigma_xe': 2.265,
                                       'B_he': 0.005, 'B_ne': 0.032, 'B_xe': 0.015, 't_total': 12}
        return look_up_table

    def generate_tipping_pulse(self):
        t0_ne = self.pulse_paras['t0_ne']
        B_he, B_ne, B_xe = self.pulse_paras['B_he'], self.pulse_paras['B_ne'], self.pulse_paras['B_xe']
        omega_he, omega_ne, omega_xe = self.omega['H'], self.omega['N'], self.omega['X']
        sigma_he = self.pulse_paras['sigma_he']
        sigma_ne = self.pulse_paras['sigma_ne']
        sigma_xe = self.pulse_paras['sigma_xe']
        gauss_A = self.gauss_A
        time_step = self.time_step
        t_total = self.pulse_paras['t_total']

        time_array1 = np.arange(0, t0_ne, time_step)
        waveform1 = np.zeros(len(time_array1))
        time_array2 = np.arange(t0_ne, t_total + time_step, time_step)
        waveform2 = np.zeros(len(time_array2))

        for i in range(len(time_array1)):
            ne = B_ne * math.sin(omega_ne * time_array1[i])
            xe = B_xe * math.sin(omega_xe * time_array1[i]) * gaussian(time_array1[i], 0, sigma_xe, gauss_A)
            he = B_he * math.sin(omega_he * time_array1[i]) * gaussian(time_array1[i], 0, sigma_he, gauss_A)
            waveform1[i] = ne + xe + he

        for i in range(len(time_array2)):
            ne = B_ne * math.sin(omega_ne * time_array2[i]) * gaussian(time_array2[i], t0_ne, sigma_ne, gauss_A)
            xe = B_xe * math.sin(omega_xe * time_array2[i]) * gaussian(time_array2[i], 0, sigma_xe, gauss_A)
            he = B_he * math.sin(omega_he * time_array2[i]) * gaussian(time_array2[i], 0, sigma_he, gauss_A)
            waveform2[i] = ne + xe + he

        time_array = np.concatenate((time_array1, time_array2), axis=None)
        pulse_B = np.concatenate((waveform1, waveform2), axis=None)
        return time_array, pulse_B

    def set_magnetic_field(self):
        By = np.zeros([len(self.time_array)])
        Bz = self.B0 * np.ones([len(self.time_array)])
        return By, Bz
