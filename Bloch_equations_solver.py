import numpy as np
import sys


class BlochEquationsSolver:
    
    def __init__(self, gamma, time_array, time_step, M0, B):
        self.gamma = gamma
        self.time_array = time_array
        self.time_step = time_step
        self.M0 = M0
        self.B = B
        self.check_length()
        Mx, My, Mz = self.runge_kutta_solver()
        self.M = [Mx, My, Mz]

    def check_length(self):
        t_l = len(self.time_array)
        bx_l = len(self.B[0])
        by_l = len(self.B[1])
        bz_l = len(self.B[2])

        if not t_l == bx_l == by_l == bz_l:
            print('The lengths of magnetic field components should equal the length of time_array')
            sys.exit()

    def be_dMx_dt(self, My, Mz, By, Bz):
        dMx_dt = self.gamma*(My*Bz - Mz*By)
        return dMx_dt
    
    def be_dMy_dt(self, Mx, Mz, Bx, Bz):
        dMy_dt = self.gamma*(Mz*Bx - Mx*Bz)
        return dMy_dt

    def be_dMz_dt(self, Mx, My, Bx, By):
        dMz_dt = self.gamma*(Mx*By - My*Bx)
        return dMz_dt

    def runge_kutta_solver(self):
        time_step = self.time_step
        Bx = self.B[0]
        By = self.B[1]
        Bz = self.B[2]
        t_l = len(self.time_array)

        Mx = np.zeros([t_l])
        My = np.zeros([t_l])
        Mz = np.zeros([t_l])

        Mx[0] = self.M0[0]
        My[0] = self.M0[1]
        Mz[0] = self.M0[2]
        for i in range(1, t_l):
            k1 = time_step*self.be_dMx_dt(My[i-1], Mz[i-1], By[i-1], Bz[i-1])
            m1 = time_step*self.be_dMy_dt(Mx[i-1], Mz[i-1], Bx[i-1], Bz[i-1])
            n1 = time_step*self.be_dMz_dt(Mx[i-1], My[i-1], Bx[i-1], By[i-1])

            k2 = time_step*self.be_dMx_dt(My[i-1]+0.5*m1, Mz[i-1]+0.5*n1, 0.5*(By[i-1]+By[i]), 0.5*(Bz[i-1]+Bz[i]))
            m2 = time_step*self.be_dMy_dt(Mx[i-1]+0.5*k1, Mz[i-1]+0.5*n1, 0.5*(Bx[i-1]+Bx[i]), 0.5*(Bz[i-1]+Bz[i]))
            n2 = time_step*self.be_dMz_dt(Mx[i-1]+0.5*k1, My[i-1]+0.5*m1, 0.5*(Bx[i-1]+Bx[i]), 0.5*(By[i-1]+By[i]))
            k3 = time_step*self.be_dMx_dt(My[i-1]+0.5*m2, Mz[i-1]+0.5*n2, 0.5*(By[i-1]+By[i]), 0.5*(Bz[i-1]+Bz[i]))
            m3 = time_step*self.be_dMy_dt(Mx[i-1]+0.5*k2, Mz[i-1]+0.5*n2, 0.5*(Bx[i-1]+Bx[i]), 0.5*(Bz[i-1]+Bz[i]))
            n3 = time_step*self.be_dMz_dt(Mx[i-1]+0.5*k2, My[i-1]+0.5*m2, 0.5*(Bx[i-1]+Bx[i]), 0.5*(By[i-1]+By[i]))

            k4 = time_step*self.be_dMx_dt(My[i-1]+m3, Mz[i-1]+n3, By[i], Bz[i])
            m4 = time_step*self.be_dMy_dt(Mx[i-1]+k3, Mz[i-1]+n3, Bx[i], Bz[i])
            n4 = time_step*self.be_dMz_dt(Mx[i-1]+k3, My[i-1]+m3, Bx[i], By[i])

            Mx[i] = Mx[i-1] + 1/6*(k1+2*k2+2*k3+k4)
            My[i] = My[i-1] + 1/6*(m1+2*m2+2*m3+m4)
            Mz[i] = Mz[i-1] + 1/6*(n1+2*n2+2*n3+n4)
        return Mx, My, Mz
