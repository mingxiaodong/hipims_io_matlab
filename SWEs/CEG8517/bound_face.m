function [zbfb,zfb,qxfb] = bound_face(bd_dir,zbfp,zfp,qxfp)
%
% Purpose: to impose proper boundary conditions.
%
% bd_dir -- location of the boundary: eastern = 1 ; western = 2
zbfb=zbfp;
if bd_dir == 1 % eastern boundary
    zfb=zfp; qxfb=qxfp;
elseif bd_dir == 2 %western boudarny
    zfb=zfp; qxfb=qxfp;
end
return
end

