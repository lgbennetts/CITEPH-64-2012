% function Param = ParamDef3d_Oceanide(GeomDisks,RIGID)
%
% LB July 2013: written for experiments

function Param = ParamDef3d_Oceanide(GeomDisks,RIGID)

if ~exist('GeomDisks','var'); GeomDisks=[0,0,0.495,33e-3]; end

% Floes are rigid:
if ~exist('RIGID','var'); RIGID=1; end

% Number of disks
Param.Np = size(GeomDisks,1);

% Accelaration due to gravity (in m\,s^{-2})
Param.g = 9.81*ones(Param.Np,1);            

% Density of fluid (in kg\,m^{-3})
Param.rho_0 = 1025*ones(Param.Np,1);           

%%%%% Choose to give draught 1.8/100 %%%%%

% Densities of disks (in kg\,m^{-3})
Param.rho = 1.8*Param.rho_0/3.3; %500*ones(Param.Np,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Poisson's ratios
Param.nu = 0.3*ones(Param.Np,1);

% Young modulii (in MPa)
if ~RIGID
 Param.E = 750*1e6*ones(Param.Np,1);
else
 Param.E = 750*1e18*ones(Param.Np,1);
end
    
% Draughts (in m)
Param.draft = GeomDisks(:,4).*Param.rho./Param.rho_0; 

% Flexural rigidities (in MPa\,m^3)
Param.D = Param.E.*GeomDisks(:,4).^3./(12*(1-Param.nu.^2));

% Scaled rigidity 
Param.beta = Param.D./(Param.g.*Param.rho_0);

% Reflection coefficient (at high frequency)
Param.Rb = 0.1;

return