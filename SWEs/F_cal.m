function [F1,F2] = F_cal(zb_L, zb_R, z_L, z_R, qx_L, qx_R)
% * purpose: to calculate the X direction flux
g = 9.81; %gravity acceleration
%****central differential solutions
zbf = (zb_L+zb_R)/2; zf = (z_L+z_R)/2; qxf = (qx_L+qx_R)/2;
hf = zf-zbf;
F1 = qxf;
F2 = qxf^2/hf+0.5*g*(zf^2-2*zf*zbf);
end