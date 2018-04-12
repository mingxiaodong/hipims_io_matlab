function [gradient] = grad(Up,Uu,Ud)
%
global dx
%
% Purpose: to gradient of flow information in a cell.
%
gradient=(Uu-Ud)/(2.0*dx);
return
%
%-------------------------------------------------------------------- %
end

