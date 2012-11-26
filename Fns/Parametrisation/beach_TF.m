function R = beach_TF(Rb,omega)

% R = beach_TF(omega)
%
% This function describes the transfer function of the beach, that is its
% reflection coefficient (in amplitudes). The TF behaves as a decaying
% exponential function of the frequency. 

A_b = 1-Rb;     % Percentage anttenuated at high frequency
c_b = 5;        % Calibration parameter

R = A_b*(exp(-c_b*omega) - 1) + 1;