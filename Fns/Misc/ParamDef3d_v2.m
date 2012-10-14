function Param = ParamDef3d_v2(Geom, thicks)

% Param = ParamDef(Geom)
%
% This function simply defines the constants of the problem. The vector
% Geom in argument contains all the scaled dimensions of geometry of the
% problem. The plate's dimensions will allow us to calculate the plate's
% non-dimensionalized parameters.

g=9.81*ones(1,Geom(2));                                 % Accelaration due to gravity
rho_water=1025*ones(1,Geom(2));                         % Density of water


rho_ice=500*ones(1,Geom(2));%922.5;                % Density of ice
nu=0.3*ones(1,Geom(2));                            % Poisson's ratio of ice
E=750*10^6*ones(1,Geom(2));%6*10^9;                % Young modulus of ice
draft=thicks.*rho_ice./rho_water;      % Draught of the plate
D=E.*(thicks*Geom(1)).^3./(12*(1-nu.^2));  % flexural rigidity
beta=D./(rho_water.*g*Geom(1)^4);                    % stiffness
gamma=draft;                                       % mass


Param=[g.' rho_water.' draft.' beta.' gamma.' rho_ice.' nu.' E.' D.'];
    