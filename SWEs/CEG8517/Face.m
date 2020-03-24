function [zbf,zf,qxf] = Face(nmod,ndir,i)
%
global M dx zb z qx zm qxm
%
% Purpose: to calculate the face values.
%
% nmod -- =1 means 1st step of RK2; = 2 means 2nd step of RK2
% ndir -- direction of the interface: eastern = 1 ; western = 2
% i -- cell i
if nmod == 1
    zp=z(i); qxp=qx(i);
else
    zp=zm(i); qxp=qxm(i);
end
zbp=zb(i);
if i+1 > M %boundary
    [zbu,zu,qxu] = bound(1,zbp,zp,qxp);
else
    if nmod == 1
        zu=z(i+1); qxu=qx(i+1);
    else
        zu=zm(i+1); qxu=qxm(i+1);
    end
    zbu=zb(i+1);
end
if i-1 < 1 %boundary
    [zbd,zd,qxd] = bound(2,zbp,zp,qxp);
else
    if nmod == 1
        zd=z(i-1); qxd=qx(i-1);
    else
        zd=zm(i-1); qxd=qxm(i-1);
    end
    zbd=zb(i-1);
end
if ndir == 1
    [grad_z] = grad(zp,zu,zd);
    zf=zp+0.5*dx*grad_z;
    [grad_qx] = grad(qxp,qxu,qxd);
    qxf=qxp+0.5*dx*grad_qx;
    zbf=(zbp+zbu)/2.0D0;
else
    [grad_z] = grad(zp,zu,zd);
    zf=zp-0.5*dx*grad_z;
    [grad_qx] = grad(qxp,qxu,qxd);
    qxf=qxp-0.5*dx*grad_qx;
    zbf=(zbp+zbd)/2.0D0;
    10
    Handout_code.m 25/04/2012
end
return
end

