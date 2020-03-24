function [E, N] = ll2os( lat, lon )
% [E, N] = LL2OS( LAT, LON ) returns UK National Grid coodinates E, N 
% (easting and northing) for points given by latitude and longitude 
% Convert latitude/longitude => OS National Grid Reference points 
% algorithm and constants adapted from D00659 v2.1 Dec 2010 of
% http://www.ordnancesurvey.co.uk/oswebsite/gps/docs/A_Guide_to_Coordinate_Systems_in_Great_Britain.pdf
% (c) Michael Fourman 2012 

  phi    = lat .* (pi/180); % convert arguments to radians
  lambda = lon .* (pi/180);
  
  a = 6377563.396;
  b = 6356256.909;              % Airy 1830 major & minor semi-axes
  F0 = 0.9996012717;            % NatGrid scale factor on central meridian
  phi0 = 49*pi/180;             % \phi_0 
  lambda0 = -2*pi/180;          % \lambda_0
  N0 = -100000;                 % NatGrid true origin 
  E0 =  400000;                 % northing & easting of true origin, metres
  e2 = 1 - b*b/(a*a);           % eccentricity squared
  n = (a-b)/(a+b);              % C1
  
  e2sin2phi = e2 * sin(phi).^2;
  nu  = a * F0 * (1-e2sin2phi).^(-0.5);
  rho = a * F0 * (1 - e2).*(1-e2sin2phi).^(-1.5);
  
  eta2 = (nu./rho) - 1 ;         %C2
  
  M = b .* F0 .* ( ...
      (1 + n + 1.25 .* (n^2 + n^3)) .* (phi - phi0) - ...
      3 .* (n + n^2 + (0.875 * n^3)) .* sin(phi - phi0) .* cos(phi + phi0) + ...
      1.875 .* (n^2 + n^3) .* sin(2 .* (phi - phi0)).* cos(2 .* (phi + phi0)) - ...
      (35/24) .* n^3 .* sin(3 .* (phi - phi0)) .* cos(3 .* (phi + phi0)) ...
      ) ;                       %C3
  I   = M + N0;
  II  = (nu / 2) .* sin(phi) .* cos(phi);
  III = (nu / 24) .* sin(phi) .* cos(phi).^3 .* (5 - tan(phi).^2 + (9 * eta2));
  IIIA = (nu./720) .* sin(phi) .* cos(phi).^5 .* ...
      (61 - 58 .* tan(phi).^2 + tan(phi).^4);
  IV = nu .* cos(phi);
  V  = (nu ./ 6) .* cos(phi).^3 .* (nu./rho - tan(phi).^2);
  VI = (nu ./ 120) .* cos(phi).^5 .* ...
      (5 - 18 .* tan(phi).^2 + tan(phi).^4 + (14 - 58 .* tan(phi).^2).* eta2);
  N = I + II.* (lambda - lambda0).^2 + IIIA .* (lambda - lambda0).^4 ; %C4
  E = E0 + IV .* (lambda - lambda0) + V .* (lambda - lambda0).^3 + VI .* (lambda - lambda0).^5;  %C5
end

