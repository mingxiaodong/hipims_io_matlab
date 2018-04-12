function [zb_b,z_b,qx_b] = bound(bd_dir,zbp,zp,qxp)
% purpose: to impose proper boundary conditions
% bd_dir: location of the boundary; eastern = 1, western = 2;
zb_b = zbp;
if bd_dir == 1 %eastern boundary
    z_b = zp; qx_b = qxp;
elseif bd_dir ==2 %western boundary
    z_b = 65; qx_b = 0;
end

end