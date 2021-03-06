# Bloch Equations Solver
This is a Bloch equation solver for Nuclear Magnetic Resonance (NMR) using Runge-Kutta methods.

### MATLAB Version

dMx.m, dMy.m and dMz.m are the three Bloch equations. Bloch_equation_solver.m solves the Bloch equations using Runger-Kutta method.

First, You need to pass the gyromagnetic ratio *gamma*. Then you need to pass the *time_array* and *time_step*. 
*Time_array* must be equal spaced and *time_step* is the spacing. You also need to pass three magnetic field arrays 
*Bx*, *By* and *Bz* which contain the magnitude of magnetic field in x, y and z directions at the time specified by 
*time_array*. Last but not least, you also need to pass the intial values of magnetization *Mx0*, *My0* and *Mz0*.

three_tipping_pulse.m is an example based on this solver. I designed a pulse to rotate Ne, Xe and He into xy plane 
simultaneously.

### Python Version
BlochEquationsSolver class in Bloch_equations_solver.py solves the Bloch equations with Runge-Kutta methods. The 
parameters you need to pass are *gamma*, *time_array*, *time_step*, *M0* and *B*. In addition to same requirements as 
MATLAB Version, *M0 = [Mx0, My0, Mz0]* and *B = [Bx, By, Bz]*. 
