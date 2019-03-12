function [d_Mz] = dMz(gamma, Mx, My, Bx, By)
% dMz/dt in Bloch equations

d_Mz = gamma*(Mx*By - My*Bx);

end

