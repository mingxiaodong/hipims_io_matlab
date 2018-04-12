function [zbfb,zfb,qxfb] = bound_face(bd_dir,zbfp,zfp,qxfp)
zbfb=zbfp;
if bd_dir==1
    zfb = zfp; qxfb = qxfp;
elseif bd_dir==2
    zfb = zfp; qxfb = qxfp;
end
end