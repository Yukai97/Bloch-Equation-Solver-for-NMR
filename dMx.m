function [d_Mx] = dMx(gamma, My, Mz, By, Bz)
% dMx/dt in Bloch equaitons

d_Mx = gamma*(My*Bz - Mz*By);

end

