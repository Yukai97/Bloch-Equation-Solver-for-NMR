clear;
clc;

%parameters
gamma_He = 2*pi*32.43409966;   %MHz/T
gamma_Ne = gamma_He/9.649;
gamma_Xe = 2*pi*11.77717;
B0 = 2*pi*15.0148/gamma_He;
omega_He = gamma_He*B0;
omega_Ne = gamma_Ne*B0;
omega_Xe = gamma_Xe*B0;

%time array and magnetic field 
t_total = 10;
time_step = 0.0002;
time_array = 0:time_step:t_total;

gaussian = @(t,t0,sigma,A)A.*exp(-(t-t0).^2/(2*sigma.^2));
Ne_BxA = 0.02;
Xe_BxA = 0.015;
He_BxA = 0.005;
gaussianA = 1;
t0_Ne = 6.3;
t0_Xe = 0;
t0_He = 0;
sigma_Ne = 0.91;
sigma_Xe = 2.26;
sigma_He = 2.468;

index_Ne = int64(t0_Ne/time_step + 1);
Ne_Bx1 = Ne_BxA*sin(omega_Ne*time_array(1:index_Ne-1));
Ne_Bx2 = Ne_BxA*sin(omega_Ne*time_array(index_Ne:end)).*gaussian(time_array(index_Ne:end),t0_Ne,sigma_Ne,gaussianA);
Ne_Bx = [Ne_Bx1, Ne_Bx2];
Xe_Bx = Xe_BxA*sin(omega_Xe*time_array).*gaussian(time_array,t0_Xe,sigma_Xe,gaussianA);
He_Bx = He_BxA*sin(omega_He*time_array).*gaussian(time_array,t0_He,sigma_He,gaussianA);
Bx = Ne_Bx + Xe_Bx + He_Bx;

By = zeros(length(time_array),1);
Bz = B0*ones(length(time_array),1);

%Initial values
Ne_Mx0 = 0;
Ne_My0 = 0;
Ne_Mz0 = 1;
Xe_Mx0 = 0;
Xe_My0 = 0;
Xe_Mz0 = 1;
He_Mx0 = 0;
He_My0 = 0;
He_Mz0 = 1;

[Ne_Mx, Ne_My, Ne_Mz] = Bloch_equation_solver(gamma_Ne, time_array, time_step, Bx, By, Bz, Ne_Mx0, Ne_My0, Ne_Mz0);
[He_Mx, He_My, He_Mz] = Bloch_equation_solver(gamma_He, time_array, time_step, Bx, By, Bz, He_Mx0, He_My0, He_Mz0);
[Xe_Mx, Xe_My, Xe_Mz] = Bloch_equation_solver(gamma_Xe, time_array, time_step, Bx, By, Bz, Xe_Mx0, Xe_My0, Xe_Mz0);
plot(time_array, Ne_Mz, time_array, Xe_Mz, time_array, He_Mz)




