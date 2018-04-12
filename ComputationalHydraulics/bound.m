function [zb_B,z_B,qx_B] = bound(bd_dir,zbp,zp,qxp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% to impose proper boundary conditions
%bd_dir -- location of the boundary 1 east; 2 west
zb_B = zbp;
if bd_dir==1
    z_B = 2;
    qx_B = qxp;
elseif bd_dir==2
    z_B = zp;
    qx_B = 4.42;
end
end

