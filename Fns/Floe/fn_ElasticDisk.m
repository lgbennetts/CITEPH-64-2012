% function fn_ElasticDisk
%
% modified from directory Arb_Floe_Pool/fn_Circ_Floe_Nov09 (Otago, 2009)
%
% calculate linear scattering of monochromatic wave by floating 
% thin-elastic disk using eigenfunction matching method
% see Peter et al (IJOPE, 2003) or Bennetts & Williams (JFM, 2010) 
%
% nb disk centred at origin
%
% nb. defined: Phi = (g/i\omega)*phi(x,y,z)*exp(-i*omega*t) velocity potential
%              xi  = eta(x,y)*exp(-i*omega*t)  surface elevation
%     gives dynamic surface condition: phi(z=0) = eta 
%
% INPUT
%
% Param = field: g (gravity), rho_0 (water density), rho (disk density)
%                nu (Poisson's ratio), E (Young's mod), draft,
%                D (flexural rigid), beta (scaled flex rigid),
%                N (vertical modes trav + ev)
% GeomDisk = [0,0,radius, thickness] 
% forceinput = {forcing type, value}
%  -> forcing type = 'p' wave period; 'k' wavenumber water; 'i' wavenumber ice
% bed = depth of fluid
% th_vec = discrete angular spectrum
% RIGID = make the rigidity effectively infinite (0=off, 1=on)
%
% OUTPUT
%
% E0 = int_{0}^{2\pi} E dtheta where ...
% E = energy radiated in unit time per unit angle (Meylan & Squire, JGR `96
%                                                  Meylan et al JGR `97)

%%
function [E0, E, E_Fourier] = fn_ElasticDisk(Param, GeomDisk, Forcing, bed, th_vec, COMM)

Tol_vec(1) = 1e-16; % - Real root error - %
Tol_vec(2) = 1e-16; % - Imag root error (init) - %
Tol_vec(3) = 1e-4;  % - Tol on const soln - %
Tol_vec(4) = 1e-1;  % - Az_Dim tol - %
Tol_vec(5) = 1e-3;  % - Energy error - %

if ~exist('COMM','var'); COMM=1; end

if ~exist('RIGID','var'); RIGID=1; end

if COMM
 path(path,'../../EXTRA_MATLAB_Fns')
end

if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.5,1e-1]; end
if ~exist('Param','var'); Param = ParamDef3d(GeomDisk,RIGID); 
    Param = ModParam_def(Param,0,0); end

draught = Param.draft; 
radius = GeomDisk(3);
Vert_Modes = Param.Ndtm;

if ~exist('bed','var'); bed=2; end

if ~exist('Forcing','var'); Forcing = Force_def(Param.g,bed,'wlength',4*radius); end

if ~exist('th_vec','var'); th_vec=linspace(-pi,pi,101); end

% DimG_vec = [DimG, Gegenbauer yes=1/no=0]
%  ->  DimG = The number of Gegenbauer Modes used (>0)   
%  ->  Gegenbauer yes/no = use the Gegenbauer polys at vertical interface (0=no)
if ~exist('DimG_vec','var'); DimG_vec=[Vert_Modes,0]; end

parameter_vector = [Param.rho_0, Param.rho, Param.nu, Param.E, Param.g];

% Normalise

Norm_fac = 1;

thickness = GeomDisk(4)/Norm_fac; draught = draught/Norm_fac; %lam = lam/Norm_fac;
radius = radius/Norm_fac; bed = bed/Norm_fac; freq = 2*pi*Forcing.f/sqrt(Norm_fac);

al = draught; % - mass
be = Param.beta; % - flex rigid

parameter_vector = [parameter_vector, ...
  (freq^2/Param.g)/Norm_fac, bed, draught, thickness, bed-draught, al, be];

%% Define matricies 

kk = zeros(Vert_Modes,1); k0 = zeros(Vert_Modes,1); % - wvnos in ice (kk) and water (k0)
wt = zeros(Vert_Modes,1); wt_0 = zeros(Vert_Modes,1); % - weight (normalising fn) attached to each vert mode

for loop_Dim=1:Vert_Modes
 kk(loop_Dim) = GetRootsMMA_PWC(parameter_vector, loop_Dim, Tol_vec);
 wt(loop_Dim) = weight_PWC(parameter_vector, kk(loop_Dim));
 k0(loop_Dim) = GetRootsMMA_FS_PWC(parameter_vector(1,:), loop_Dim, Tol_vec);
 wt_0(loop_Dim) = weight_0_PWC(bed, k0(loop_Dim));   
end

[mu_0, mu_1] = mu_new_PWC(parameter_vector, Vert_Modes, kk, wt); % - apx complex roots (Dimension dependent)

% s_Modes = The limt of the Fourier expansion in horizontal plane (any nat num)
%           total number = 2*s_Modes+1

Az_Dim_vec = 2*Tol_vec(4)*ones(3,1);
Az_Dim = 0;
count = 0;
while max(Az_Dim_vec)>Tol_vec(4)    
    Az_Dim_vec(count+1) = abs(besselj(Az_Dim,k0(1)*radius));
    Az_Dim=Az_Dim+1; count = count+1; count = mod(count,length(Az_Dim_vec));
end

clear Az_Dim_vec count

Az_Dim=Az_Dim-2; s_Modes = Az_Dim; clear Az_Dim 

% DimG = DimG_vec(1);
% 
% mat_A = matrix_A_PWC(parameter_vector, Vert_Modes, kk, wt, parameter_vector(10)); % - in ice 
% mat_A0 = matrix_A_PWC(parameter_vector, Vert_Modes, k0, wt_0, parameter_vector(7)); % - in water

[mat_A0,mat_A,mat_W0, mat_W] = fn_JumpMats(parameter_vector,...
    Vert_Modes, DimG_vec(1), k0, kk, wt_0, wt,DimG_vec(2));

mat_A0 = diag(mat_A0);

C_mat = mat_C_dr_PWC(parameter_vector, length(kk), kk, mu_0, mu_1,...
     wt, parameter_vector(12)); % - matrix that controls motion in ice

mat_W0_T = conj(mat_W0'); mat_W_T = conj(mat_W');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve via EMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Psi = zeros(Vert_Modes, 2*s_Modes+1); Phi = zeros(Vert_Modes, 2*s_Modes+1); 

% -- Version 1 -- %

Avec = zeros(Vert_Modes, 2*s_Modes+1); Bvec = zeros(Vert_Modes, 2*s_Modes+1); 

v_01 = C_mat(1:Vert_Modes, Vert_Modes+1:Vert_Modes+2); %clear C_mat

% -- Version 1.1 -- %

for loop_Az=-s_Modes:s_Modes
    
%% - Incident wave

Inc = (1i^loop_Az)*exp(abs(imag(k0(1)*radius))).*eye(Vert_Modes,1);

%% Bend & Shear
% NB. scale besselj by exp(-abs(imag(z))) & besselh by exp(-iz)

[bend_mu0, shear_mu0] = ...
    BendShear3dAxi_fn( parameter_vector, [1, loop_Az, 1], mu_0*radius );
[bend_mu1, shear_mu1] = ...
    BendShear3dAxi_fn( parameter_vector, [1, loop_Az, 1], mu_1*radius );
[bend_K, shear_K] = ...
    BendShear3dAxi_fn( parameter_vector, [1, loop_Az, 1], diag([kk*radius]) );

% Bessel fns

Amp_mu(1,1:Vert_Modes) = -C_mat(Vert_Modes+1, 1:Vert_Modes)*...
    (bend_mu1*shear_K - shear_mu1*bend_K);
Amp_mu(2,1:Vert_Modes) = C_mat(Vert_Modes+1, 1:Vert_Modes)*...
    (bend_mu0*shear_K - shear_mu0*bend_K);

Amp_mu = -Amp_mu./(bend_mu0*shear_mu1 - bend_mu1*shear_mu0);
    
Amp_mu_Mat(:,:,loop_Az+s_Modes+1) = diag(exp(-imag([mu_0 mu_1])*radius))*Amp_mu;
   
%% Components of full solution: Bessel & Hankel fns 

J_ice = besselj(loop_Az, kk.*radius, 1);
Jz_ice = Bessel_dz(@besselj, loop_Az, kk.*radius, 1);

Jmu_ice = besselj(loop_Az, [ mu_0 mu_1 ].*radius, 1);
Jzmu_ice = Bessel_dz(@besselj, loop_Az, [ mu_0 mu_1 ].*radius, 1);

J_water = besselj(loop_Az, k0.*radius, 1);
Jz_water = Bessel_dz(@besselj, loop_Az, k0.*radius, 1);

H_water = besselh(loop_Az, 1, k0.*radius, 1);
Hz_water = Bessel_dz(@besselh, loop_Az, k0.*radius, 1, 1);

%% Set up soln in water:

Phi_In = diag(J_water);
Phiz_In = diag(k0.*Jz_water);

Phi_Out = diag(H_water);
Phiz_Out = diag(k0.*Hz_water);

%% Set up soln from ice:

% - fn psi = %

Psi_S = diag(J_ice) + v_01*diag(Jmu_ice)*Amp_mu;

Psiz_S = diag(kk.*Jz_ice) + v_01*diag([mu_0, mu_1].*Jzmu_ice)*Amp_mu;

w1 = 1:Vert_Modes; w2 = Vert_Modes+1:Vert_Modes+2;

% w_S = C_mat(Vert_Modes+1,w1).*J_ice.' +...
%     C_mat(Vert_Modes+1,w2)*diag(Jmu_ice)*Amp_mu;
% 
% w_In = C_mat(Vert_Modes+1,w1).*J_ice.' +...
%     C_mat(Vert_Modes+1,w2)*diag(Hmu_ice)*Amp_mu_J;
% w_Out = C_mat(Vert_Modes+1,w1).*H_ice.' +...
%     C_mat(Vert_Modes+1,w2)*diag(Hmu_ice)*Amp_mu_H;
% wht_In = C_mat(Vert_Modes+2,w1).*J_ice.' +...
%     C_mat(Vert_Modes+2,w2)*diag(Hmu_ice)*Amp_mu_J;
% wht_Out = C_mat(Vert_Modes+2,w1).*H_ice.' +...
%     C_mat(Vert_Modes+2,w2)*diag(Hmu_ice)*Amp_mu_H;
% 
% wz_In = C_mat(Vert_Modes+1,w1).*(kk.*Jz_ice).' +...
%     C_mat(Vert_Modes+1,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_J;
% wz_Out = C_mat(Vert_Modes+1,w1).*(kk.*Hz_ice).' +...
%     C_mat(Vert_Modes+1,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_H;
% wzht_In = C_mat(Vert_Modes+2,w1).*(kk.*Jz_ice).' +...
%     C_mat(Vert_Modes+2,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_J;
% wzht_Out = C_mat(Vert_Modes+2,w1).*(kk.*Hz_ice).' +...
%     C_mat(Vert_Modes+2,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_H;

%% Solve the system

% - following are Dim x DimG:

MatB_u0 = (Phiz_Out\(diag(mat_A0)\mat_W0)); 
MatB_u = (Psiz_S\(mat_A\mat_W)); 

% - following is Dim x Dim:

MatB_I0 = -Phiz_Out\Phiz_In;

% - following are DimG x Dim:

MatA_I0 = mat_W0_T*(Phi_In + Phi_Out*MatB_I0);

% - following are DimG x DimG:

MatA_u0 = mat_W0_T*Phi_Out*MatB_u0;
MatA_u = mat_W_T*Psi_S*MatB_u;

% - find u in terms of the inc amps (DimG x Dim)

MatU_I0 = (MatA_u - MatA_u0)\MatA_I0;

u_vec = MatU_I0*Inc;

% - convert into scattered amps - %
% - in the ice (floe) - %
Avec(:, s_Modes+loop_Az+1) = MatB_u*u_vec; 
% - in the (surrounding) water - %
Bvec(:, s_Modes+loop_Az+1) = MatB_u0*u_vec + MatB_I0*Inc;

% - Put into vel pots and disp - %

Phi(:, s_Modes+loop_Az+1) = Phi_In*Inc + Phi_Out*Bvec(:, s_Modes+loop_Az+1);
Psi(:, s_Modes+loop_Az+1) = Psi_S*Avec(:, s_Modes+loop_Az+1);

% w(:, s_Modes+loop_Az+1) = w_S*Avec(:, s_Modes+loop_Az+1);
% wn(:, s_Modes+loop_Az+1) = wz_In*Inc + wz_Out*Bvec(:, s_Modes+loop_Az+1);
% 
% what(:, s_Modes+loop_Az+1) = wht_In*Inc + wht_Out*Bvec(:, s_Modes+loop_Az+1);
% wnhat(:, s_Modes+loop_Az+1) = wzht_In*Inc + wzht_Out*Bvec(:, s_Modes+loop_Az+1);

end

%% Energy check

% In far-field:
% \phi=\sqrt{2/pi/k0/r}exp(i*k0*r)\sum_{n=-N}^{N}Bvec(1,n)*exp(i*n*(theta-\pi/2))
%     =\sqrt{2/pi/k0/r}exp(i*k0*r)*A(theta)
% where 
% A(theta) = \sum_{n=-N}^{N}Bvec(1,n)*exp(i*n*(theta-\pi/2))

% Amplitudes of angular components of prop wave:
Bvec0 = Bvec(1,:);
% Eliminate in-built scaling of Hankel functions (see: help besselh)
Bvec0 = exp(-1i*k0(1)*radius)*Bvec0;
% Useful coefficients:
Ip = exp(1i*[-s_Modes:s_Modes]*pi/2); 

% 1. E0 = int_{0}^{2\pi}|A(theta)|^2 dtheta 
%    using orthogonality = ... 
E0 = sum(abs(Bvec0.*conj(Ip)).^2);

% 2. A(in_angle) = A(0) = \sum_{n=-N}^{N}Bvec(1,n)*exp(n*\pi/2i))
%                       = \sum_{n=-N}^{N}Bvec(1,n)*Ip(n)
% => Real{A(0)} = E1 = ...
E1 = real(sum(conj(Ip).*Bvec0));

% int_{0}^{2\pi}|A(theta)|^2 dtheta = -Real{A(inc_angle)}
if abs(E0+E1) > Tol_vec(5)
 cprintf('red',['energy error = ', num2str(E0+E1), '\n'])
end

%% Natural modes of vibration

Jmax=5;

% Displacement

res = 101;
displ=zeros(2*s_Modes+1,res);
r_vec = linspace(0,radius,res);

for loop_r=1:length(r_vec)
 count=1;
 for loop_Az=-s_Modes:s_Modes
  JJ = besselj(loop_Az, kk*r_vec(loop_r));
  JJ_mu = besselj(loop_Az, [mu_0;mu_1]*r_vec(loop_r));
    
  dum = C_mat(Vert_Modes+1,w1)*diag(JJ) + ...
   C_mat(Vert_Modes+1,w2(1))*JJ_mu(1)*Amp_mu_Mat(1,:,s_Modes+1+loop_Az) + ...
   C_mat(Vert_Modes+1,w2(2))*JJ_mu(2)*Amp_mu_Mat(2,:,s_Modes+1+loop_Az);

  displ(count,loop_r) = dum*Avec(:,s_Modes+1+loop_Az);
  count = count+1;
 end
end
clear count dum

if 1
 th=0;
 
 evec = exp(1i*[-s_Modes:s_Modes]*th);
 evec = reshape(evec,1,2*s_Modes+1);
 full_displ = zeros(1,length(r_vec));
 
 for loop_r=1:length(r_vec)
  full_displ(loop_r) = evec*displ(:,loop_r);
 end
 
 full_displ = full_displ/cosh(k0(1)*bed)/wt_0(1);
 
 figure; plot(r_vec,real(full_displ),'r'); hold on
 plot(r_vec,real(exp(1i*k0(1)*r_vec)),'ko');
 clear evec full_displ
end
     
% Inner-products

A_even=zeros(s_Modes+1,Jmax+1);
A_odd =zeros(s_Modes+1,Jmax+1);

%[freq,shape] = Circ_plate_nat_freq(nu,N,J)

for loop_Az=0:s_Modes
 for loop_j=0:Jmax   
  [wnj,Nnj] = fn_NatModes(loop_Az,loop_j,Param.nu,r_vec,radius);
  if loop_Az~=0
   A_even(loop_Az+1,loop_j+1) = ...
       pi*fn_Trap(r_vec.*(displ(s_Modes+1+loop_Az,:)+...
             displ(s_Modes+1-loop_Az,:)),wnj,r_vec)/Nnj;
   A_odd(loop_Az+1,loop_j+1) = ...
       pi*fn_Trap(r_vec.*(displ(s_Modes+1+loop_Az,:)-...
             displ(s_Modes+1-loop_Az,:)),wnj,r_vec)/Nnj;
  else
   A_even(loop_Az+1,loop_j+1) = ...
       2*pi*fn_Trap(r_vec.*displ(s_Modes+1,:),wnj,r_vec)/Nnj;
  end
 end
end 

% Ensure even solution
if max(max(abs(A_odd)))>1e-5
 disp('A_odd = ')
 disp(abs(A_odd))
end

% Ensure rigid plate (if specified)
if RIGID
 dum_A=A_even; dum_A(1,1)=0; dum_A(2,1)=0;  
 if max(max(abs(dum_A)))>1e-3
  disp('A_even = ')
  disp(abs(A_even))
 end
 clear dum_A
end

%% Output

% FF_Amp = Far-field amplitude (function of theta) 
% E ~ |A|^2 -> scattered energy (function of energy)
% E0 = int_{0}^{2\pi} E dtheta

% nb. incamp not required as axisymm problem
% incang = 0;
% Normalise inc wave energy to unity with incamp = sqrt(2/rho0*g)
% incamp = sqrt(2/Param.rho_0*Param.g);
% but this simply cancels with below, so leave !!!!
% (CC Mei book Chap 1)

% A(theta) = \sum_{n=-N}^{N}Bvec(1,n)*exp(i*n*(theta-\pi/2)) as above
FF_Amp = (Bvec0.*conj(Ip))*exp(1i*[-s_Modes:s_Modes].'*th_vec);
% Set the inc wave amplitude to 1 at this point:
FF_Amp = FF_Amp/wt_0(1)/cosh(k0(1)*bed);

E = 2*(abs(FF_Amp).^2)/k0(1)/pi;

E0 = E0/2/pi;

% Fourier basis representation
% Array E(i,j) such that E = sum_{p.q=-N}^{N}E(p,q)*exp(i*(p-q)*theta)

Bvec0 = (Bvec0.*conj(Ip));

E_Fourier = zeros(2*s_Modes+1,2*s_Modes+1);

for loop_s1=1:2*s_Modes+1    
 for loop_s2=1:2*s_Modes+1
  E_Fourier(loop_s1,loop_s2)=Bvec0(loop_s1)*conj(Bvec0(loop_s2));
 end
end

E_Fourier = 2*E_Fourier./k0(1)./pi;
 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          %%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mat_A0,mat_A,mat_W0,mat_W] = ...
    fn_JumpMats(pv,Vert_Dim,DimG,k0,kk,wt_0,wt,METH)

%% -- coeffs of the JCs at r=R -- %

mat_A = matrix_A_PWC(pv, Vert_Dim, kk, wt, pv(10));
mat_A0 = matrix_A_PWC(pv, Vert_Dim, k0, wt_0, pv(7));

if METH==0
 %%% - Old way (vertical modes) - %%%
 %[mat_W0, mat_W] = jump_W_PWC(pv, Vert_Dim, DimG, kk, k0, kk, wt_0, wt);
 [mat_W0, mat_W] = jump_W_PWC_wtd(pv, Vert_Dim, kk, k0, kk, wt, wt_0, wt);
else
 %%% - New way (Gegenbauer) - %%%
 if pv(9)==0
    [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, pv(7), pv(8), 0.5, DimG);
 else
    [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, pv(7), pv(8), 1/6, DimG);
 end
end

%% -- Put jumps together in big matrices - %

% Big_AW0 = zeros((2*Az_Dim+1)*Vert_Dim, (2*Az_Dim+1)*DimG); Big_AW = Big_AW0;
% Big_VT0 = Big_AW0.'; Big_VT = Big_VT0;
% 
% v1 = 1:Vert_Dim;
% v2 = 1:DimG;
% 
% for loop_Az = 1:2*Az_Dim+1
%     Big_AW0(v1, v2) = mat_A0\mat_W0;
%     Big_AW(v1, v2) = mat_A\mat_W;
%     
%     Big_VT0(v2,v1) = mat_W0.';
%     Big_VT(v2,v1) = mat_W.';
%     
%     v1 = v1+Vert_Dim; v2 = v2+DimG;
% end

return

% NATURAL MODES 
%
% see Montiel 2012 (Thesis, p.150)
%
% wnj(r,th) = ( Jn(Knj*r/R)+Cnj*In(Knj*r/R) )*cos(n*th) (n>0)
% wnj(r,th) = ( Jn(Knj*r/R)+Cnj*In(Knj*r/R) )           (n=0)
% wnj(r,th) = ( Jn(Knj*r/R)+Cnj*In(Knj*r/R) )*sin(n*th) (n<0)
%
% Nnj = normalisation
% nb. Knj=K-nj, Cnj=(-1)^n*C-nj and Nnj=N-nj
% asymptotic behaviour:
% Knj = (pi/2)(n+2j)   (Leissa, 1969)
% Thus use interval [(pi/2)(n+2(j-1)), (pi/2)(n+2j)] in fzero.  

function [wnj,Nnj] = fn_NatModes(n,j,nu,r_vec,R)

if and(n==0,j==0)
    
 wnj = 1+0*r_vec;
 
elseif and(abs(n)==1,j==0)
    
 wnj = r_vec;
 
else
    
 params(1) = n; params(2)=nu;
 
 int = [(pi/2)*(n+2*(j-1)), (pi/2)*(n+2*j)];

 Knj = fn_bisection(@fn_disprelNM, int, params);

%  fj = (Knj.^3).*Bessel_dz(@besselj,n,Knj) + (n^2)*(1-nu)*(...
%     Knj.*Bessel_dz(@besselj,n,Knj) - besselj(n,Knj) );
% 
%  fi = (Knj.^3).*Bessel_dz(@besseli,n,Knj) - (n^2)*(1-nu)*(...
%     Knj.*Bessel_dz(@besseli,n,Knj) - besseli(n,Knj) );
% 
%  Cnj = fj/fi; clear fi fj

 mj = (Knj.^2).*besselj(n,Knj) + (1-nu)*(...
    Knj.*Bessel_dz(@besselj,n,Knj) - (n^2)*besselj(n,Knj) );

 mi = (Knj.^2).*besseli(n,Knj) - (1-nu)*(...
    Knj.*Bessel_dz(@besseli,n,Knj) - (n^2)*besseli(n,Knj) );

 Cnj = mj/mi; clear mi mj

 wnj = besselj(n,Knj*r_vec/R) + Cnj*besseli(n,Knj*r_vec/R);
 
end

if and(n==0,j==0)
 Nnj=pi*(R^2);
elseif and(abs(n)==1,j==0)
 Nnj = pi*(R^4)/4;
else
 Nnj = besselj(n,Knj)^2 - besselj(n-1,Knj)*besselj(n+1,Knj) + ...
     (Cnj/Knj)*(besselj(n,Knj)*besseli(n-1,Knj) - ...
      besselj(n-1,Knj)*besseli(n,Knj) ) + ...
     (Cnj^2)*(besseli(n,Knj)^2 - besseli(n-1,Knj)*besseli(n+1,Knj) );
 Nnj = pi*Nnj*(R^2)/2;   
end

return

% 

function D=fn_disprelNM(K,params)

n = params(1); nu=params(2);

fj = (K.^3).*Bessel_dz(@besselj,n,K) + (n^2)*(1-nu)*(...
    K.*Bessel_dz(@besselj,n,K) - besselj(n,K) );

fi = (K.^3).*Bessel_dz(@besseli,n,K) - (n^2)*(1-nu)*(...
    K.*Bessel_dz(@besseli,n,K) - besseli(n,K) );

mj = (K.^2).*besselj(n,K) + (1-nu)*(...
    K.*Bessel_dz(@besselj,n,K) - (n^2)*besselj(n,K) );

mi = (K.^2).*besseli(n,K) - (1-nu)*(...
    K.*Bessel_dz(@besseli,n,K) - (n^2)*besseli(n,K) );

for loop=1:length(K)
 D(loop)=det([fj(loop),mj(loop);fi(loop),mi(loop)]);
end

return

% Fabien's bisection method

function out = fn_bisection(fn_name,int,params)

guess_min=int(1); guess_max=int(2); clear int

if guess_min == 0
    guess_min = guess_min+0.01;
end

ans1 = guess_min; ans2 = guess_max;

while abs(ans2 - ans1) > 1e-9
    fval = feval(fn_name,(ans1+ans2)/2, params);
    if fval*feval(fn_name,ans1,params) < 0
        ans2 = (ans1+ans2)/2;
    else
        ans1 = (ans1+ans2)/2;
    end
end

out = ans2;

return

% 

function I = fn_Trap(f,g,x)

dx = x(2)-x(1);

f = reshape(f,1,length(x));
g = reshape(g,1,length(x));

I = f.*g;

I(1) = I(1)/2; I(end) = I(end)/2;

I = dx*sum(I);

return

%% OLD OUTPUT CODE

% - Displacement - % 
%
% ang=0;
% res = 50; w_vec=zeros(1,2*res-1);
% dum_w=zeros(1,res);
% r_vec = linspace(0,radius,res);
% 
% for loop_Az=-s_Modes:s_Modes
%     
%     JJ = besselj(loop_Az, kk*r_vec);
%     JJ_mu = besselj(loop_Az, [mu_0;mu_1]*r_vec);
%     
%     dum_vec = C_mat(Vert_Modes+1,w1)*JJ +...
%         C_mat(Vert_Modes+1,w2)*(diag(Amp_mu_Mat(:,s_Modes+1+loop_Az))*Jmu_ice);
% 
%     dum_w = dum_w + ...
%         exp(1i*loop_Az*(ang+pi))*w(:,s_Modes+loop_Az+1)*dum_vec;
% 
%     w_vec(res:2*res-1) = w_vec(res:2*res-1) + ...
%         exp(1i*loop_Az*ang)*w(:,s_Modes+loop_Az+1)*dum_vec;
% end
% 
% loop_r = 1;
% for loop_rr=res:-1:2
%     w_vec(loop_r) = dum_w(loop_rr); loop_r=loop_r+1;
% end
% 
% r_vec = linspace(-radius,radius,2*res-1);
%
% % - Edge of floe - %
%  
% Avec0 = Psi(1,:);
% EI = sum(abs(Avec0).^2);
% EI = EI/2/pi/radius;
%  
% % - Other stuff - %
% 
% k_out(1) = k0(1); k_out(2) = kk(1);
% Avec1 = Avec(1,:);

% % - Kochin function - %
% 
% cprintf('red','Dont use -> not debugged -> use Far-field amplitude instead\n')
%
% fac = Param.rho_0*(freq^3)/8/pi/Param.g;
% 
% for loop_th=1:length(th_vec)
%  E(loop_th) = fn_Kochin(pi+th_vec(loop_th),Psi,Amp_mu,k0(1),kk,...
%      [mu_0,mu_1],wt,bed-draught,radius);
% end
% 
% E = fac*(abs(E).^2);

%% BELOW IS NOT DEBUGGED - ALSO NOT REQUIRED

% % Kochin function
% 
% function H = fn_Kochin(tau,amps,amps_mu,k0,kk,mus,wt,dep,rad)
% 
% H = 0;
% 
% NM = size(amps);
% 
% for loopN = 1:NM(1) % vertical modes
%  kap = kk(loopN); wt_n = wt(loopN);
%  M = -(NM(2)-1)/2;
%  a_mu = amps_mu(:,loopN);
%  for loopM = 1:NM(2) % azimuthal modes
%   H = H + amps(loopN,loopM)*wt_n*(k0*cosh(kap*dep) - kap*sinh(kap*dep))*...
%       (fn_IntBess(abs(M),kap,k0,rad) + ...
%       a_mu(1)*fn_IntBess(abs(M),mus(1),k0,rad) + ...
%       a_mu(2)*fn_IntBess(abs(M),mus(2),k0,rad) ) ...
%                                             *exp(1i*M*(tau-pi/2));
%   M = M + 1;
%  end
% end
%      
% H = 2*pi*H;
% 
% return
% 
% % int_{0}^{R}[besselj(n,a*r)*besselj(-n,b*r)*r]dr
% % using wolframalpha.com
% 
% function I = fn_IntBess(n,a,b,R) 
% 
% if a~=b
% 
%  if n~=0
%  
%   r = R;
%  
%   I = r*( b*besselj(n,a*r)*besselj(n-1,b*r) - ...
%       a*besselj(n,b*r)*besselj(n-1,a*r) );
%   
%   r = 0;
%  
%   I = I - r*( b*besselj(n,a*r)*besselj(n-1,b*r) - ...
%       a*besselj(n,b*r)*besselj(n-1,a*r) );
%  
%   I = ((-1)^n)*I/(a^2 - b^2);
%  
%  else
%     
%   r = R;
%  
%   I = -r*( b*besselj(0,a*r)*besselj(1,b*r) - ...
%       a*besselj(0,b*r)*besselj(1,a*r) );
%   
%   r = 0;
%  
%   I = I + r*( b*besselj(0,a*r)*besselj(1,b*r) - ...
%       a*besselj(0,b*r)*besselj(1,a*r) );
%  
%   I = ((-1)^n)*I/(a^2 - b^2);
%  
%  end
%  
% else
%     
%  if n~=0
%  
%   r = R;
%  
%   I = 0.5*(r^2)*( besselj(n,a*r)^2 - ...
%       besselj(n+1,a*r)*besselj(n-1,a*r) );
%   
%   r = 0;
%  
%   I = I - 0.5*(r^2)*( besselj(n,a*r)^2 - ...
%       besselj(n+1,a*r)*besselj(n-1,a*r) );
%  
%   I = ((-1)^n)*I;
%  
%  else
%     
%   r = R;
%  
%   I = 0.5*(r^2)*( besselj(0,a*r)^2 + besselj(1,b*r)^2 );
%   
%   r = 0;
%  
%   I = I - 0.5*(r^2)*( besselj(0,a*r)^2 + besselj(1,b*r)^2 );
%  
%   I = ((-1)^n)*I;
%  
%  end  
%  
% end
%  
% return

