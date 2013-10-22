% function Bvec = fn_CircScat(s_Modes, forceinput, radius) 
%
% DESCRIPTION: Calculates Diffraction Transfer Matrix for circular scatterer 
% with Helmholtz equation in exterior & sound-hard bcs, i.e. phi_n=0 on r=R
%
% INPUTS:
%
% s_Modes = The limt of the Fourier expansion in horizontal plane (any nat num)
%           total number = 2*s_Modes+1
% forceinput = {forcing type, value}
%  -> forcing type = 'l' wavelength; 'k' wavenumber
% radius = radius of floe (>0)
%
% L Bennetts Sept 2013 / Adealide 
%
% REVISION HISTORY:
%
% - copied & modified from Directional_spectrum/Main_Files/fn_CircScat.m
% - modified from Arb_Floe_Pool/fn_Circ_Floe_Nov09.m 
  
function Bvec = fn_CircScat(s_Modes, forceinput, radius) 
    
if ~exist('s_Modes','var'); s_Modes = 2; end
if ~exist('forceinput','var'); forceinput = {'l',1}; end
if ~exist('radius','var'); radius = 0.5; end

%% - Define geom parameters - %

if forceinput{1}=='l'
    k0 = 2*pi/forceinput{2};
elseif forceinput{1}=='k'
    k0 = forceinput{2};
end

Vert_Modes = 1;

%% - Solve - %
Phi = zeros(Vert_Modes, 2*s_Modes+1); 

%% -- Version 1 -- %

Bvec = zeros(2*s_Modes+1); 
       
%% -- Bessel & Hankel fns -- %

%J_water = besselj(loop_Az, k0.*radius, 1);
Jz_water = Bessel_dz(@besselj, -s_Modes:s_Modes, k0.*radius, 1);

% - Note that the scaling of besselh is exp(-i*kr)

%H_water = besselh(loop_Az, 1, k0.*radius, 1);
Hz_water = Bessel_dz(@besselh, -s_Modes:s_Modes, k0.*radius, 1, 1);

%% -- Set up soln in water -- %

%Phi_In = diag(J_water);
Phiz_In = diag(k0.*Jz_water);

%Phi_Out = diag(H_water);
Phiz_Out = diag(k0.*Hz_water);

%% -- Solve -- %

Bvec = -Phiz_Out\Phiz_In;

% - Put into vel pots and disp - %

%Phi = Phi_In*Inc + Phi_Out*Bvec(:, s_Modes+loop_Az+1);

%% - Energy relation
 
Inc = ( (1i.^[-s_Modes:s_Modes])*exp(abs(imag(k0(1)*radius))) ).';

Bvec0 = (Bvec*Inc).';
% - NOTE: the following removes the scaling of besselh
Bvec0 = exp(-1i*k0(1)*radius)*Bvec0;
% Bvec0 = Bvec0.*exp([-s_Modes:s_Modes]*pi/2i);
Ip = exp(1i*[-s_Modes:s_Modes]*pi/2); %inc_amp*wt(1)*

E0 = sum(abs(Bvec0).^2);
E1 = real(sum(conj(Ip).*Bvec0));
if abs(E0+E1) > 1e-3
    display(['energy error', num2str(E0+E1)])
end

% E0 = E0/2/pi;
% 
% H_amps = Bvec0.';
% 
% %% - Directional Spectrum
% 
% theta=pi/3;
% 
% [dum_Sp(1,:), Se, T, w] = besselh_dirspec(-s_Modes,theta);
% Sp = Bvec0(1)*dum_Sp;
% 
% for n=-s_Modes+1:s_Modes
%  [dum_Sp(n+s_Modes+1,:), Se, T, w] = besselh_dirspec(n,theta);
% end
% 
% Sp = Bvec0*dum_Sp;
% 
% %% - Far-field amplitude
%  
% incang = pi/6;
%  
% th_vec = linspace(-pi,pi,101);
%  
% Bvec0 = Bvec0.*exp([-s_Modes:s_Modes]*pi/2i);
%  
% FF_Amp = Bvec0*exp(1i*[-s_Modes:s_Modes].'*(th_vec-incang));
  
return
