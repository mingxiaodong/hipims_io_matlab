function [F1,F2] = F_cal(zb_L,z_L,qx_L,zb_R,z_R,qx_R)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% to calculate the x-direction flux
g = 9.81;
zbf = (zb_L+zb_R)/2;
zf = (z_L+z_R)/2;
qxf = (qx_L+qx_R)/2;
hf = zf-zbf;
F1 = qxf;
F2 = qxf^2/hf +0.5*g*(zf^2-2.0*zf*zbf);
end

