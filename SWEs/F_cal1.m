function [F1,F2] = F_cal1(zbf,z_L, z_R, qx_L, qx_R)
% * purpose: to calculate the X direction flux
g = 9.81; %gravity acceleration
%**to calculate h and u, and a = sqrt(g*h)
h_L = z_L-zbf; u_L = qx_L/h_L; a_L = sqrt(g*h_L); 
h_R = z_R-zbf; u_R = qx_R/h_R; a_R = sqrt(g*h_R);

%**to calculate wave speed SL and SR
u_star = (u_L+u_R)/2 + a_L-a_R;
h_star = (1/g)*((a_L+a_R)/2 + (u_L-u_R)/4)^2;
s_L = min(u_L-a_L, u_star-sqrt(g*h_star)); 
s_R = min(u_R+a_R, u_star+sqrt(g*h_star));

% FL and FR are two-element vectors
F_L1 = qx_L; F_L2 = u_L*qx_L+0.5*g*z_L*(z_L-2.0*zbf);
F_R1 = qx_R; F_R2 = u_R*qx_R+0.5*g*z_R*(z_R-2.0*zbf);

if s_L>=0.0
    F1 = F_L1; F2 = F_L2;
else
    if s_R>=0
        F1=(s_R*F_L1 - s_L*F_R1 + s_L*s_R*( z_R- z_L) ) / (s_R-s_L);
        F2=(s_R*F_L2 - s_L*F_R2 + s_L*s_R*(qx_R-qx_L) ) / (s_R-s_L);
    else
        F1 = F_R1; F2 = F_R2;
    end
end

end