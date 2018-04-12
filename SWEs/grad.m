function [gradient] = grad(Up,Uu,Ud)
global dx
gradient=(Uu-Ud)/(2.0*dx);
end