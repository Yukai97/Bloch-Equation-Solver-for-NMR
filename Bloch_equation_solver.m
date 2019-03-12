function [Mx, My, Mz] = Bloch_equation_solver(gamma, time_array, time_step, Bx, By, Bz, Mx0, My0, Mz0)
%Solving Bloch equation solver woth Runge-Kutta method

L1 = length(time_array);
L2 = length(Bx);
L3 = length(By);
L4 = length(Bz);

if L1 ~= L2
    error('The length of Bx should equal the length of time_array');
elseif L1 ~= L3
    error('The length of By should equal the length of time_array');
elseif L1 ~= L4
    error('The length of Bz should equal the length of time_array');
end

Mx = zeros(1,L1);
My = zeros(1,L1);
Mz = zeros(1,L1);

Mx(1) = Mx0;
My(1) = My0;
Mz(1) = Mz0;

%Runge-Kutta method
for i = 2:L1
    k1 = time_step*dMx(gamma,My(i-1),Mz(i-1),By(i-1),Bz(i-1));
    m1 = time_step*dMy(gamma,Mx(i-1),Mz(i-1),Bx(i-1),Bz(i-1));
    n1 = time_step*dMz(gamma,Mx(i-1),My(i-1),Bx(i-1),By(i-1));
    
    k2 = time_step*dMx(gamma,My(i-1)+0.5*m1,Mz(i-1)+0.5*n1,0.5*(By(i-1)+By(i)),0.5*(Bz(i-1)+Bz(i)));
    m2 = time_step*dMy(gamma,Mx(i-1)+0.5*k1,Mz(i-1)+0.5*n1,0.5*(Bx(i-1)+Bx(i)),0.5*(Bz(i-1)+Bz(i)));
    n2 = time_step*dMz(gamma,Mx(i-1)+0.5*k1,My(i-1)+0.5*m1,0.5*(Bx(i-1)+Bx(i)),0.5*(By(i-1)+By(i)));
    
    k3 = time_step*dMx(gamma,My(i-1)+0.5*m2,Mz(i-1)+0.5*n2,0.5*(By(i-1)+By(i)),0.5*(Bz(i-1)+Bz(i)));
    m3 = time_step*dMy(gamma,Mx(i-1)+0.5*k2,Mz(i-1)+0.5*n2,0.5*(Bx(i-1)+Bx(i)),0.5*(Bz(i-1)+Bz(i)));
    n3 = time_step*dMz(gamma,Mx(i-1)+0.5*k2,My(i-1)+0.5*m2,0.5*(Bx(i-1)+Bx(i)),0.5*(By(i-1)+By(i)));
    
    k4 = time_step*dMx(gamma,My(i-1)+m3,Mz(i-1)+n3,By(i),Bz(i));
    m4 = time_step*dMy(gamma,Mx(i-1)+k3,Mz(i-1)+n3,Bx(i),Bz(i));
    n4 = time_step*dMz(gamma,Mx(i-1)+k3,My(i-1)+m3,Bx(i),By(i));
    
    Mx(i) = Mx(i-1) + 1/6*(k1+2*k2+2*k3+k4);
    My(i) = My(i-1) + 1/6*(m1+2*m2+2*m3+m4);
    Mz(i) = Mz(i-1) + 1/6*(n1+2*n2+2*n3+n4);
end
    

end

