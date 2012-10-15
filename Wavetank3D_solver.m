function [Coef_1P,Coef_1M,Coef_3P,Coef_3M] = Wavetank3D_solver(...
    Af,k,data_disks,x,W,H,omega,N,M,Rw,Rb)

% [Coef_1P,Coef_1M,Coef_3P,Coef_3M] = Wavetank3D_solver(...
%     Af,k,data_disks,x,W,H,omega,N,M,Rw,Rb)
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
% data_disks contains all the data required to solve the problem of 
% multiple elastic disks in a channel (to be defined by Luke).
% x = [x0 x1 x2 x3] as defined on the schematic above.
% W is the wavetank width.
% H is the depth.
% omega is the radian frequency.
% N+1 is the number of vertical modes (free surface region).
% M+1 is the number of transversal modes (tank modes).
% Rw is the reflection coefficient of the wavemaker.
% Rb is the reflection coefficient of the beach.


%% Definition of the modal weighting
Y_weight = 1;
Z_weight = Vert_Weight(k,H);


%% Wavemaker scattering at x0: Coef_1P = Sw * Coef_1M + F.
[Sw,F] = Wavemaker3D_solver(Af,k,Y_weight,Z_weight,H,omega,M,N,Rw,W,...
    x(1),x(2));


%% Scattering by a group of elastic disks in a channel located in [x1 x2]
% [Coef_1M Coef_2P].' = Sd * [Coef_1P Coef_2M].'
Sd = DisksInChannel_solver(k,data_disks,Y_weight,Z_weight,x(2),x(3),W,H,...
    omega,N,M);


%% Beach scattering at x3: Coef_2M = Sb * Coef_2P
Sb=Beach3D_solver(k,M,N,omega,x(end),x(end-1));


%% Solution of the wave tank problem

% Mapping matrices
Map1P_to_3P = (eye((N+1)*(M+1)) - ...
    Sd((N+1)*(M+1)+1:end,(N+1)*(M+1)+1:end)*Sb)\...
    Sd((N+1)*(M+1)+1:end,1:(N+1)*(M+1));

Map1P_to_1M = Sd(1:(N+1)*(M+1),1:(N+1)*(M+1)) + ...
    Sd(1:(N+1)*(M+1),(N+1)*(M+1)+1:end)*Sb*Map1P_to_3P;

MapF_to_1P=eye((N+1)*(M+1))-Sw*Map1P_to_1M;

% Vectors of amplitudes
Coef_1P = MapF_to_1P\F;
Coef_1M = Map1P_to_1M*Coef_1P;
Coef_3P = Map1P_to_3P*Coef_1P;
Coef_3M = Sb*Coef_3P;