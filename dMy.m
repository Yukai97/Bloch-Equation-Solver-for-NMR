function [d_My] = dMy(gamma, Mx, Mz, Bx, Bz)
% dMy/dt in Bloch equations

d_My = gamma*(Mz*Bx - Mx*Bz);

end

