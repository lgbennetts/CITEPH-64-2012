function [Coef_1P, Coef_1M, Coef_3P, Coef_3M] = Wavetank3D_solver(...
    Af, k, Z_weight, v_vec, Rm, Tm, Rp, Tp, W, H, x, omega, N, M, P, ...
    Rb, Rw)

% [Coef_1P,Coef_1M,Coef_3P,Coef_3M] = Wavetank3D_solver(...
%     Af,k,Z_weight,v_vec,Rm,Tm,Rp,Tp,W,H,x,omega,N,M,P,Rw)
%
% Cross-sectional view of the wave tank:
%
%  wavemaker                                                     artificial
%        /                                                 z=0      beach
%       /----------------------|                   |-------------------|-
%      /   -->         <--     | region with disks |    -->      <--   |
%z=Zh /  Phi_1P       Phi_1M                          Phi_3P   Phi_3M  |
%     |                        |                   |                   |
%     |                                                                |
%     |        (1)             |        (2)        |        (3)        |
%     |                                                                |
%     |__________z=-H__________|___________________|___________________|_
%    x=x0                     x=x1                x=x2                x=x3
%
% The output vectors are the vectors of potential amplitudes for the waves
% travelling as specified in the above schematic. 
%
% Af is the amplitude of the wavemaker.
% k is the vector of wavenumbers in free surface region.
% Z_weight contains weights of vertical modes. 
% v_vec contains the horizontal wavenumbers scaled by k.
% [Rm, Tp; Tm, Rp] is the scattering matrix of disks region.
% W is the wavetank width.
% H is the depth.
% x = [x0 x1 x2 x3] as defined on the schematic above.
% omega is the radian frequency.
% N+1 is the number of vertical modes (free surface region).
% M+1 is the number of transversal modes (tank modes).
% P is the number of propagating horizontal modes.
% Rw is the reflection coefficient of the wavemaker.

%% Longitudinal wavenumbers
alpha = zeros((N+1)*(M+1),1);
for m = 0:M
    alpha(m*(N+1)+1:(m+1)*(N+1)) = k.'.*v_vec(:,m+1);
end


%% Wavemaker scattering at x0: Coef_1P = Sw * Coef_1M + F.
[Sw,F] = Wavemaker3D_solver(Af, k, Z_weight, alpha, 1, H, omega, M, N, ...
    Rw, W, x(1), x(2));


%% Beach scattering at x3: Coef_2M = Sb * Coef_2P
Sb = Beach3D_solver(alpha, M, P, N, Rb, omega, x(end), x(end-1));


%% Solution of the wave tank problem

% Mapping matrices
Map1P_to_3P = (eye((N+1)*(M+1)) - Rp*Sb)\Tm;

Map1P_to_1M = Rm + Tp*Sb*Map1P_to_3P;

MapF_to_1P=eye((N+1)*(M+1))-Sw*Map1P_to_1M;

% Vectors of amplitudes
Coef_1P = MapF_to_1P\F;
Coef_1M = Map1P_to_1M*Coef_1P;
Coef_3P = Map1P_to_3P*Coef_1P;
Coef_3M = Sb*Coef_3P;