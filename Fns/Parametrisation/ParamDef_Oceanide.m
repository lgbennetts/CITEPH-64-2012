% function Param = ParamDef_Oceanide(GeomDisks,RIGID)
%
% LB July 2013: written for experiments

function Param = ParamDef_Oceanide(RIGID,Np)

%%Param is structure with attributes and default values:
%%          Np: 1                         no of disks 
%%   thickness: 0.033000000000000         thickness [m]
%%           g: 9.810000000000000         gravity [m s^{-2}]
%%       rho_0: 1025                      water density [kg m^{-3}]
%%         rho: 5.590909090909091e+02     floe density [kg m^{-3}]
%%          nu: 0.300000000000000         Poisson's ratio
%%           E: 1.000000000000000e+10     Youngs modulus [Pa]
%%       draft: 0.018000000000000         draft, rho/rho_0*thickness [Pa m^3]
%%           D: 3.290934065934066e+04     flexural rigidity
%%        beta: 3.272851561059214         scaled rigidity, D/rho/g [m^4]
%%         bed: 3.100000000000000         water depth [m]
%%  MIZ_length: 5                         length of MIZ [m]
%%   floe_diam: 0.990000000000000         floe diameter [m]

% Floes are rigid:
if ~exist('RIGID','var'); RIGID=10; end

% Number of disks
if ~exist('Np','var')
 Param.Np = 1;
else
 Param.Np = Np;
end

% thickness
Param.thickness = 33e-3;

% Accelaration due to gravity (in m\,s^{-2})
Param.g = 9.81*ones(Param.Np,1);            

% Density of fluid (in kg\,m^{-3})
Param.rho_0 = 1025*ones(Param.Np,1);           

%%%%% Choose to give draught 1.8/100 %%%%%

% Densities of disks (in kg\,m^{-3})
Param.rho = 1.8*Param.rho_0/3.3; %500*ones(Param.Np,1);

% Poisson's ratios
Param.nu = 0.3*ones(Param.Np,1);

% Young's modulii (in Pa)
if ~RIGID
 Param.E = 10e10*ones(Param.Np,1);
else
 Param.E = RIGID*(10^9)*ones(Param.Np,1);
end
    
% Draughts (in m)
Param.draft = Param.thickness.*Param.rho./Param.rho_0; 

% Flexural rigidities (in MPa\,m^3)
Param.D = Param.E.*Param.thickness.^3./(12*(1-Param.nu.^2));

% Scaled rigidity 
Param.beta = Param.D./(Param.g.*Param.rho_0);

% fluid depth (open/equilibrium)

Param.bed = 3.1;

% MIZ length [m]

Param.MIZ_length = 5;

% diameter of floes [m]

Param.floe_diam = .99;

return
