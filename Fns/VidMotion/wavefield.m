% [field] = wavefield(prop,val,h)
%
% LJ YIEW
% Created on  Feb 2014
% Last edited Oct 2016
%
% Calculates the corresponding wave properties of a wave field for a given
% input parameter, according to the dispersion relation.
%
% INPUTS:
% prop  = 'lambda' - wavelength [m]
%         'f'      - frequency [Hz]
%         'omega'  - angular frequency [rad/s]
%         'T'      - wave period [s]
%         'k'      - wavenumber
%         'sigma'  - omega^2/g 
% val   = enter value
% h     = water depth [m]
%
% OUTPUTS:
% [field] = 6x2 cell e.g.
%           'Frequency'            [0.1103]
%           'Angular Frequency'    [0.6929]
%           'Wave Period'          [9.0675]
%           'Wavelength'           [    20]
%           'Wavenumber'           [0.3142]
%           'Water Depth'          [0.5000]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [field] = wavefield(prop,val,h)

g = 9.81;

if prop == 'lambda'
 
 lambda = val;
 k      = 2*pi/lambda;
 omega  = sqrt(g*k*tanh(k*h));
 f      = omega/(2*pi);
 T      = 1/f;
 sigma  = omega^2/g;
 
elseif prop == 'k'
 
 k      = val;
 lambda = 2*pi/k;
 omega  = sqrt(g*k*tanh(k*h));
 f      = omega/(2*pi);
 T      = 1/f;
 sigma  = omega^2/g;
  
elseif prop == 'f'
 
 f      = val;
 omega  = 2*pi*f;
 sigma  = omega^2/g;
 T      = 1/f;
 k      = abs(fzero(@(k) k*tanh(k*h)-sigma,0));
 lambda = 2*pi/k;
 
elseif prop == 'T'
 
 T      = val;
 f      = 1/T;
 omega  = 2*pi*f;
 sigma  = omega^2/g;
 k      = abs(fzero(@(k) k*tanh(k*h)-sigma,0));
 lambda = 2*pi/k;
  
elseif prop == 'omega'
 
 omega  = val;
 f      = omega/(2*pi);
 T      = 1/f;
 sigma  = omega^2/g;
 k      = abs(fzero(@(k) k*tanh(k*h)-sigma,0));
 lambda = 2*pi/k;
 
elseif prop == 'sigma'
 
 sigma  = val;
 omega  = sqrt(sigma*g);
 f      = omega/(2*pi);
 T      = 1/f;
 k      = abs(fzero(@(k) k*tanh(k*h)-sigma,0));
 lambda = 2*pi/k;
 
end

field = [cellstr('Frequency')         f;
         cellstr('Angular Frequency') omega;
         cellstr('Wave Period')       T;
         cellstr('Wavelength')        lambda;
         cellstr('Wavenumber')        k;
         cellstr('Water Depth')       h];
 

 

 
 
 
 
 
