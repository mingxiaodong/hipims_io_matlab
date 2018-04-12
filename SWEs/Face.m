function [zbf,zf,qxf] = Face(nmod,ndir,i)
%Purpose: to calculate face values
global M dx zb z qx zm qxm

if nmod == 1
    zp = z(i); qxp = qx(i);
else
    zp =zm(i); qxp = qxm(i);
end
zbp = zb(i);

if i+1 > M %boundary
    [zbu,zu,qxu] = bound(1,zbp,zp,qxp);
else
    if nmod == 1
        zu = z(i+1); qxu = qx(i+1);
    else
        zu = zm(i+1); qxu = qxm(i+1);
    end
    zbu = zb(i+1);
end

if i-1 < 1
    [zbd,zd,qxd] = bound(2,zbp,zp,qxp);
else
    if nmod ==1
        zd = z(i-1);qxd = qx(i-1);
    else
        zd=zm(i-1); qxd=qxm(i-1);
    end
    zbd=zb(i-1);
end

if ndir == 1
    
    [grad_z] = grad(zp,zu,zd);
    zf = zp+0.5*dx*grad_z;
    [grad_qx] = grad(qxp,qxu,qxd);
    qxf=qxp+0.5*dx*grad_qx;
    zbf=(zbp+zbu)/2.0;
else
    [grad_z]=grad(zp,zu,zd);
    zf=zp-0.5*dx*grad_z;
    [grad_qx]=grad(qxp,qxu,qxd);
    qxf=qxp-0.5*dx*grad_qx;
    zbf=(zbp+zbd)/2.0;
end
end
