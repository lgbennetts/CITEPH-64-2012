% function fn_2dFloe(TEST,fortyp,lam0,conc,Ens_size,COMM)
%
% DESCRIPTION: THIS FUNCTION WILL CALCULATE REFLECTION & TRANSMISSION
% COEFFICIENTS FOR FLOES OF UNIFORM THICKNESS: THE LENGTH OF THE FLOES
% IS EITHER (A) LONG; OR (B) RANDOMLY SELECTED FROM A UNIFORM DIST.
%
% nb. Attn model used in the WIM as at July 2013 (+ floe length dep)
% nb. vertical modes cosh(k(z+h))/cosh(kh) fixed
%
% INPUTS:
%
% TEST = which problem
%        option 1 = 'Oceanide'    
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam_vec = a vector of forcing fortyp
% conc = concentration of ice cover
% Ens_size = how many members of ensemble
%            only necessary when long floe limit off
% COMM = comments on or off
% Nd   = number of vertical modes
%
% OUTPUTS:
%
% TT    = Transmitted energy
%
% EXTRA VARIABLES:
%
% parameter_vector = [liq_dens, ice_dens, poiss, youngs, grav, freq^2/g]
% visc_rp          = viscosity parameter (default = 0)

function out = ...
 fn_2dFloe(fortyp,lam0,Param,outputs,LONG,COMM,Ens_size)

%% Inputs & prelims

if ~exist('DO_PLOT','var'); DO_PLOT=0; end

if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('lam0','var');   lam0=1/2.5; end

if ~exist('LONG','var'); LONG=0; end
if ~exist('COMM','var'); COMM=1; end

if ~exist('RIGID','var'); RIGID=10; end

if ~exist('Param','var'); Param = ParamDef_Oceanide(5); 
    Param = ModParam_def(Param,1,5,0,0); end

if ~LONG; if ~exist('Ens_size','var'); Ens_size = 1; end; end

if ~exist('outputs','var'); outputs='transmitted energy heave pitch'; end

depth = Param.bed;

Forcing = Force_def(Param.g(1), depth, fortyp, lam0);

Vert_Dim = Param.Ndtm;

thickness = Param.thickness;

floe_length = Param.floe_diam;

DimG = Vert_Dim; clear Dims

draught = Param.draft; al = draught;

be = Param.beta;

fq = (2*pi*Forcing.f)^2/Param.g;

parameter_vector = [Param.rho_0, Param.rho, Param.nu, Param.E, ...
    Param.g, fq]; 

%visc_rp = 0; %1e-1; %1e-1; % 

%% FREE-SURF ROOTS

%Roots0 = roots_rp0(parameter_vector, 0, Vert_Dim-1, 0, depth, 0);
%Roots0=Roots0.';

Roots0 = zeros(Vert_Dim,1); 

Roots0(1)   = fn_RealRoot([depth*parameter_vector(6)], ...
          'fn_ReDispRel_water', 'fn_UppLimReal_water', 1e-16)/depth; 

for loop_Dim = 2:Vert_Dim
 Roots0(loop_Dim) = 1i*fn_ImagRoot_water(loop_Dim-1, depth, ...
  parameter_vector(6), 1e-16); 
end

mat_A0 = wtd_cosh_cosh(Vert_Dim, Roots0, depth);     

%% Begin calculations:

if COMM
  cprintf([0.3,0.3,0.3],'------------------------------------------\n')
  cprintf([0.3,0.3,0.3],'----------   START: 2d Attn    -----------\n')
end

if COMM
 cprintf([0.3,0.3,0.3],['>> ' fortyp ' = ' num2str(lam0) '\n'])
 if ~LONG 
  cprintf([0.3,0.3,0.3],['>> floe length = ' num2str(floe_length) '\n'])
  cprintf([0.3,0.3,0.3],['>> ensemble = ' num2str(Ens_size) '\n'])
 else
  cprintf([0.3,0.3,0.3],'>> long floe limit\n')
 end
 cprintf([0.3,0.3,0.3],['>> ' num2str(thickness) ' thick\n'])
 cprintf([0.3,0.3,0.3],['>> rigidity = ' sprintf('%0.5g',Param.E) '\n'])
 cprintf([0.3,0.3,0.3],['>>> ' int2str(Vert_Dim) ' vertical modes\n'])
 cprintf([0.3,0.3,0.3],['>>> lam0/k0 = ' num2str(2*pi/Roots0(1)) '/' num2str(Roots0(1)) '\n'])
end
     
clear Param
     
%%% ONE INTERFACE 

[Rm0,Tm0,Rp0,Tp0,Roots] = ...
            fn_WaterIce(parameter_vector, Vert_Dim, DimG,... 
            Roots0, mat_A0, thickness, depth, draught, al, be, 0);
           
if COMM; cprintf([0.3,0.3,0.3],['>>> lam/k   = ' ...
  num2str(2*pi/Roots(1)) '/' num2str(Roots(1)) '\n']); end           
        
if LONG %%% LONG FLOE LIMIT %%%
 
 r11 = Rm0(1,1); 
 %alpha = -2*log(1-abs(r11)^2);
 
 TT = (1-abs(r11)^2)^2; RR = [];

else %%% NO LONG FLOE LIMIT %%% 
            
 sd = min(floe_length,pi/Roots(1));         

 for loop=1:Ens_size   
    
  if Ens_size==1
   dum_fl = floe_length;
  else
   dum_fl = floe_length + 2*sd*(rand - 0.5);
  end
        
  [Rm,Tm,Rp,Tp] = ...
            fn_IndFloe(dum_fl, Rm0, Tm0, Rp0, Tp0, Roots);
        
  r_vec(loop) = abs(Rm(1,1))^2;
  t_vec(loop) = abs(Tm(1,1))^2;      
        
 end
 
 Rm = mean(r_vec); Tm = mean(t_vec);

 TT = Tm; RR = Rm;

end

%% PLOT:

if or(DO_PLOT,or(~isempty(strfind(outputs,'heave')),...
  ~isempty(strfind(outputs,'pitch'))))
 
 x_res = 501;
 
 [wts,xx,fac] = fn_NumInt(x_res,-floe_length/2,floe_length/2);
 
 xx_ext = (floe_length/2+2*pi/Roots0(1))*linspace(-1,1,x_res);
 
 v1 = 1:Vert_Dim+2; v2 = v1 + Vert_Dim+2;
 
 am = eye(Vert_Dim,1); bm = zeros(Vert_Dim,1);
 
 %%% Amplitudes in floe
 
 El = exp(1i*Roots*floe_length);
 
 dum_M(v1,v1)=eye(Vert_Dim+2); dum_M(v2,v2)=eye(Vert_Dim+2);
 dum_M(v1,v2)=-Rp0*diag(El); dum_M(v2,v1)=-Rp0*diag(El);
 
 dum_v(v1,1) = Tm0*am;       dum_v(v2,1) =Tm0*bm;
 
 dum_amps = dum_M\dum_v;
 
 %%% The floe profile
 
 eta=zeros(x_res,1);
 for loop_x=1:length(xx)
  eta(loop_x) = ((Roots.*tanh(Roots*(depth-draught))).')*...
   (diag(exp(1i*Roots*(xx(loop_x)+floe_length/2)))*dum_amps(v1) ...
   + diag(exp(-1i*Roots*(xx(loop_x)-floe_length/2)))*dum_amps(v2));
 end
 
 eta = eta/fq;
 
 %%% Plot
 
 if DO_PLOT
  figure(fn_getfig)
  
  h1 = subplot(2,1,1); hold on; h2 = subplot(2,1,2); hold on
  
  plot(h1,xx_ext,real(exp(1i*Roots0(1)*(xx_ext+floe_length/2))),'k:')
  plot(h1,xx,real(eta),'r')
  plot(h2,xx_ext,imag(exp(1i*Roots0(1)*(xx_ext+floe_length/2))),'k:')
  plot(h2,xx,imag(eta),'r')
 end
 
end

 
%% OUTPUTS & FINISH 

out_str = ' ''dummy'' '; out_val = ' 0 ';

if strfind(outputs,'reflected energy')
 out_str = [out_str '; ''reflected energy'' '];
 out_val = [out_val '; RR'];
end

if strfind(outputs,'transmitted energy')
 out_str = [out_str '; ''transmitted energy'' '];
 out_val = [out_val '; TT'];
end

if strfind(outputs,'heave')
 out_str = [out_str '; ''heave'' '];
 out_val = [out_val '; abs((fac/floe_length)*wts*eta)'];
end

if strfind(outputs,'pitch')
 out_str = [out_str '; ''pitch'' '];
 out_val = [out_val '; abs((12*fac/floe_length^3)*wts*(eta.*xx))'];
end

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});']) 
out(1)=[];

if COMM
 disp(['>>>> reflected energy  : ' num2str(RR)])
 disp(['>>>> transmitted energy: ' num2str(TT)])  
 cprintf([0.3,0.3,0.3],'-----------    END: 2d Attn   ------------\n')
 cprintf([0.3,0.3,0.3],'------------------------------------------\n')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUBFUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           funct to solve for individual floes                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rm,Tm,Rp,Tp] = ...
    fn_IndFloe(l_vec,Rm0,Tm0,Rp0,Tp0,kk)

No_IceInts = 1; % - aprox the roughness length
 
l_vec = [0;reshape(l_vec,No_IceInts,1);0];

%%%%%%%%%%%%%%%%%%

list_Rm{1} = Rm0;
list_Tm{1} = Tm0;
list_Rp{1} = Rp0;
list_Tp{1} = Tp0;

% - USE SYMMETRY!!!!! - %

list_Rp{No_IceInts+1} = list_Rm{1};
list_Tp{No_IceInts+1} = list_Tm{1};
list_Rm{No_IceInts+1} = list_Rp{1};
list_Tm{No_IceInts+1} = list_Tp{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   ENERGY METHOD  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Inter_Dim0 = length(list_Tm{1}(1,:));
til_Inter_Dim = length(list_Tm{1}(:,1));

dum_Mat = zeros(2*til_Inter_Dim); 
dum_Vec = zeros(2*til_Inter_Dim,Inter_Dim0+til_Inter_Dim);

v1=1:til_Inter_Dim; v2=v1+til_Inter_Dim;
u1=1:Inter_Dim0; u2=Inter_Dim0+[1:til_Inter_Dim];

w1=1:Inter_Dim0; w2=w1(end)+[1:til_Inter_Dim];
w3=w2(end)+[1:til_Inter_Dim]; w4=w3(end)+[1:til_Inter_Dim];

R0m = list_Rm{1}; T0p = list_Tp{1};
T0m = list_Tm{1}; R0p = list_Rp{1};
     
II = eye(til_Inter_Dim);

u2=u1+Inter_Dim0;

Ep = diag(exp(1i*kk*l_vec(No_IceInts+1)));

Rm = list_Rm{No_IceInts+1}; Tp = list_Tp{No_IceInts+1};
Tm = list_Tm{No_IceInts+1}; Rp = list_Rp{No_IceInts+1};

% - Version 2 - %
   
dum_Mat(v1,v1) = II; dum_Mat(v1,v2) = -R0p*Ep;
dum_Mat(v2,v2) = II; dum_Mat(v2,v1) = -Rm*Ep;

%display(abs(det(dum_Mat)))
%plot(real(det(dum_Mat)),imag(det(dum_Mat)),'bo')

dum_Vec(v1,u1) = T0m; dum_Vec(v2,u2) = Tp;

dum_Vec = dum_Mat\dum_Vec; %display([l_vec(No_IceInts+1),abs(det(dum_Mat))])

dum_R0m = R0m + T0p*Ep*dum_Vec(v2,u1);
dum_T0p = T0p*Ep*dum_Vec(v2,u2);
dum_R0p = Rp + Tm*Ep*dum_Vec(v1,u2);
dum_T0m = Tm*Ep*dum_Vec(v1,u1);

% - Version 3 - %

% v1 = 1:Inter_Dim0; v2 = v1(end) + [1:til_Inter_Dim];
% v3 = v2(end) + [1:til_Inter_Dim]; v4 = v3(end) + [1:Inter_Dim0];
% 
% II0 = eye(Inter_Dim0);
% 
% dum_Mat(v1,v1) = II0; dum_Mat(v1,v3) = -T0p*Ep;
% dum_Mat(v2,v2) = II; dum_Mat(v2,v3) = -R0p*Ep;
% dum_Mat(v3,v3) = II; dum_Mat(v3,v2) = -Rm*Ep;
% dum_Mat(v4,v4) = II0; dum_Mat(v4,v2) = -Tm*Ep;
% 
% dum_Vec(v1,u1) = R0m;
% dum_Vec(v2,u1) = T0m;
% dum_Vec(v3,u2) = Tp;
% dum_Vec(v4,u2) = Rp;
% 
% dum_Vec = dum_Mat\dum_Vec; display([l_vec(No_IceInts+1),abs(det(dum_Mat))])
% 
% dum_R0m = dum_Vec(v1,u1);
% dum_T0m = dum_Vec(v4,u1);
% dum_T0p = dum_Vec(v1,u2);
% dum_R0p = dum_Vec(v4,u2);

% -------------- %

Rm = dum_R0m; Tp = dum_T0p; Tm = dum_T0m; Rp = dum_R0p;
 
if 1
 S=[Rm(1,1),Tp(1,1);Tm(1,1),Rp(1,1)];
 if abs(abs(det(S))-1)>1e-3
   disp('error in ice Scattering matrix!!!!!!!!!!!!!!!!!')
   disp(['----->',num2str(abs(abs(det(S))-1))])
 end
end
 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          funct to solve for individual ice edge              %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rm,Tm,Rp,Tp,Roots_mat] = ...
    fn_WaterIce(parameter_vector, Vert_Dim, DimG,...
    Roots0, mat_A0, capD_vec, depth, d_vec, al_vec, be_vec, visc_rp)
    
No_Sects = 1;

% capD_vec = Di;

% d_vec = (parameter_vector(2)/parameter_vector(1))*capD_vec; %(5/6)*capD_vec;

% al_vec = (parameter_vector(2)/parameter_vector(1))*capD_vec;
% be_vec = parameter_vector(4)*(capD_vec.^3)./...
%     (12*parameter_vector(1)*(1-(parameter_vector(3)^2))*parameter_vector(5));

% - INITIALISE - %

Mu_vec = zeros(2,No_Sects);
C_mats = zeros(Vert_Dim+2,Vert_Dim+2,No_Sects);
A_mats = zeros(Vert_Dim,Vert_Dim,No_Sects);

% - CALC THE ICE COVERED ROOTS & WEIGHTS - %

% [Roots_mat, errs]  = roots_rp(parameter_vector, visc_rp, Vert_Dim-1, ...
%     capD_vec, depth - d_vec, 0);
% Roots_mat=Roots_mat.';
% if Vert_Dim>1
%  dum_perm=[1,4:Vert_Dim+2,2:3];
%  Roots_mat = Roots_mat(dum_perm,:);
% end
% if sum(sum(errs))>0, disp('There are errors in the roots.'), end

Roots_mat = zeros(Vert_Dim+2,No_Sects);
 
for loop_s=1:No_Sects
 
 pv = [parameter_vector(1:6), depth, d_vec(loop_s), capD_vec(loop_s),...
     depth-d_vec(loop_s), al_vec(loop_s), be_vec(loop_s)];
    
 Roots_mat(1,loop_s) = fn_RealRoot([1-pv(6)*al_vec(loop_s), ...
  be_vec(loop_s)/((depth-d_vec(loop_s))^4),...
  (depth-d_vec(loop_s))*pv(6)],...
  'fn_ReDispRel_ice', 'fn_UppLimReal_ice', 1e-16)/(depth-d_vec(loop_s)); 
 
 for loop_Dim = 2:Vert_Dim
  Roots_mat(loop_Dim,loop_s) = 1i*fn_ImagRoot_ice(loop_Dim-1, ...
   al_vec(loop_s), be_vec(loop_s), depth-d_vec(loop_s), ...
   pv(6), 1e-16*[1,1]);
 end

 A_mats(:,:,loop_s) = wtd_cosh_cosh(Vert_Dim, ...
     Roots_mat(1:Vert_Dim,loop_s), depth-d_vec(loop_s));
  
 [Mu_vec(1,loop_s), Mu_vec(2,loop_s)] = mu_new_PWC_visc_wtd(pv, Vert_Dim, ...
     Roots_mat(1:Vert_Dim,loop_s), A_mats(:,:,loop_s), visc_rp);
 
 C_mats(:,:,loop_s) = mat_C_dr_PWC_wtd(pv, Vert_Dim, Roots_mat(1:Vert_Dim,loop_s), ...
     Mu_vec(1,loop_s), Mu_vec(2,loop_s), A_mats(:,:,loop_s));

 clear pv
 
end

Roots_mat(Vert_Dim+[1:2],:) = Mu_vec;

Inter_Dim=Vert_Dim; cmpx_mkr=1; Inter_Dim0=Vert_Dim;

if cmpx_mkr == 1 % - Keep the cmplx rts 
    
 Roots_Inter = Roots_mat;
 Roots_Inter(Inter_Dim+1:Vert_Dim,:)=[];
 
 c_vecs = C_mats(Vert_Dim+1,[1:Inter_Dim,Vert_Dim+1,Vert_Dim+2],:);
 
 til_Inter_Dim = Inter_Dim+2;
 
else % - Eliminate the cmplx rts
    
 Roots_Inter = Roots_mat;
 Roots_Inter(Inter_Dim+1:Vert_Dim+2,:)=[];
 
 c_vecs = C_mats(Vert_Dim+1,[1:Inter_Dim],:);
 
 til_Inter_Dim = Inter_Dim;
 
end
 
Geom0 = [capD_vec(1) d_vec(1) depth al_vec(1) be_vec(1)];

[Rm, Tm, Rp, Tp] = CalcWaterIceTransition_wtd(parameter_vector, Vert_Dim, ...
    DimG, Inter_Dim0, Inter_Dim, ...
    Roots0, Roots_mat(1:Vert_Dim,1), Roots_mat(end-1:end,1),...
    C_mats(:,:,1), mat_A0, A_mats(:,:,1), Geom0, cmpx_mkr);

if and(1,visc_rp==0)
 S=[Rm(1,1),Tp(1,1);Tm(1,1),Rp(1,1)];
 if abs(abs(det(S))-1)>1e-3
   disp('error in Scattering matrix 2!!!!!!!!!!!!!!!!!')
   disp(['----->',num2str(abs(abs(det(S))-1))])
 end
end

return  

%

function parameter_vector = Param_vec(freq)

liquid_density = 1025; %1000; %  %in kg/m^3                                                                                                                                                                                                                                                                            
ice_density = 922.5; %900; %  %in kg/m^3 %                                                                                                                                                                                                                                                                             
poissons_ratio = 0.3; %1/3; %                                                                                                                                                                                         
youngs_mod = 6e9; %5000000000; 

gravity = 9.81; %1; % in m/sec^2
kappa = (freq^2)/gravity;

parameter_vector = [liquid_density, ice_density, poissons_ratio, ...
    youngs_mod, gravity, kappa];
    
return

%

function cosh_cosh = wtd_cosh_cosh(Dim, root, water_depth )

%Written by Luke
%Abbreviated by Gareth
%Rewriiten by Luke to inc a weighting

%Weighting chosen to be m_i=sech(k_i H) (06.05.08)

%Evaluates integrals of the form
%int_{-h}^{-d} cosh( root_i(it)*(z+h) ) cosh( root_j(jt)*(z+h) ) dz
%where z=-h is the seafloor, and -d is the underside of the ice

cosh_cosh(1:Dim,1:Dim,1:length(water_depth)) = complex(0);

for kt = 1:length(water_depth)
  for it = 1:Dim
    for jt = 1:Dim
      if root(it,kt) ~= root(jt,kt)
        tanh_i = tanh(root(it,kt)*water_depth(kt));
        tanh_j = tanh(root(jt,kt)*water_depth(kt));
        quot = root(it,kt)^2 - root(jt,kt)^2;
        quot = 1 / quot;
        cosh_cosh(it,jt,kt) = quot*(root(it,kt)*tanh_i - root(jt,kt)*tanh_j);
      end
      if root(it,kt) == root(jt,kt)
        if root(it,kt) == 0
          cosh_cosh(it,jt,kt) = water_depth(kt);
        else
          arg = root(it,kt)*water_depth(kt);
          tanh_i = tanh(arg);
          cosh_cosh(it,jt,kt) = (tanh_i + arg*(sech(arg)^2))/(2*root(it,kt));
        end
      end
    end
  end
 end

return

% ROOTS

function [kn errs] = roots_rp0(parameter_vector, visc_rp, Dim, h_vec, fluid_depth, check_errs_flag)

%Returns the dimensional roots of the Robinson-Palmer dispersion relation.
%Any values for period, visc_rp and h_vec may be used.

tol2 = 1e-5; %used in finding the shifted evanescent roots:critical, so don't change it 
tol = 1e-9; 
tinysteps_1 = 15;
%  tinysteps_2 = 15;
E=parameter_vector(4);%Pa
g=parameter_vector(5);%m/s^2
rho_w=parameter_vector(1);%kg/m^3
rho_ice=parameter_vector(2);%kg/m^3
nu=parameter_vector(3);
%  depth = 1000; %Used in the finite depth case
%  varpic = 5*(1/4)^.8 *real( (-1).^-.2 );


lh = length(h_vec);
kn=complex(zeros([lh, Dim+1]));
errs(1:lh,1:6)=complex(zeros([lh,6]));
  
    
for ht = 1:lh %For each thickness in the vector h_vec
  k = complex(zeros([1,Dim+1]));
  k_evan = k(1:Dim);
  om = sqrt(g*parameter_vector(6)); %2*pi/period;
  D = E* h_vec(ht)^3 / (12*(1-nu^2));
  C0 = -rho_w * om^2;
  C5 = D ;
%    L = (D/(rho_w*om^2))^0.2;
  C1r = rho_w * g  - om^2*rho_ice*h_vec(ht);
  C1i = 1i*visc_rp*om/tinysteps_1;
  C1 = C1r - 1i*visc_rp*om;
%    varpi = g / ( om^2 * L) - rho_ice * h_vec(ht) / ( rho_w*L); %Tim's non-dim param, but its useful here

  %starting values from the infinite depth case, with viscosity
  inf_rts = roots([C5 0 0 0 C1 C0]);
  [Y,I] = sort(2*pi*(angle(inf_rts)/pi <0 ) + angle(inf_rts));
  inf_rts = inf_rts(I);    %the roots, sorted by angle
  %inf_rts(3) = -inf_rts(5);
  
  %For thickness in [0.001m 35m], period in [5s 35s] and for gamma in [0,5e3] the 
  %primary root has a smaller imaginary part than the shifted, 1st quadrant root.
%   if imag(inf_rts(1)) > imag(inf_rts(2)) %Then they two roots need to be swapped
%     inf_rts(1:2) = inf_rts(2:1);
%   end

%    if abs(varpi-varpic)<0.01
%      %The dkdg method for distinguishing roo0ts doesn't always work.
%      %For varpi close to varpic use a tinysteps approach with a Newton-Raphson
%      
%      C1i = 1i*visc_rp*om/tinysteps_2;
%      
%      %starting values from the infinite depth case, zero viscosity
%      inf_rts = roots([C5 0 0 0 C1r C0]);
%      [Y,I] = sort(2*pi*(angle(inf_rts)/pi <0 ) + angle(inf_rts));
%      inf_rts = inf_rts(I);    %the roots, sorted by angle
%      inf_rts(3) = -inf_rts(5);
%  
%      for nt = 1:3;
%        for tt = 1:tinysteps_2
%          x = inf_rts(nt);
%          for ct = 1:100
%            CC1=C1r+C1i*tt;
%            %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%            f = C0 + CC1*x + C5*x^5;
%            df = CC1 + 5*C5*x^4;
%            %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%            if abs(f)<tol, break, end
%            x = x - f / df;
%          end
%          inf_rts(nt) = x;
%        end
%      end
%    else
%      %starting values from the infinite depth case, with viscosity
%      inf_rts = roots([C5 0 0 0 C1 C0]);
%      [Y,I] = sort(2*pi*(angle(inf_rts)/pi <0 ) + angle(inf_rts));
%      inf_rts = inf_rts(I);    %the roots, sorted by angle
%      inf_rts(3) = -inf_rts(5);
%     
%      
%      %It seems the derivative of the roots with respect viscosity provides a means for checking 
%      %which of the first quadrant roots is the primary and which is the complex, provided visc
%      dkdg = 1i*inf_rts(1:2).^2*L^2 ./ (4*L^5*inf_rts(1:2).^5 + 1);
%      %sort the values in inf_rts.
%      %for low viscosities angle alone distinguishes the two roots (and the dk/dg test doesn't work). 
%      if visc_rp*om > 0.1  %For larger viscosities use the dk/dg test
%        if varpi < varpic
%          %dkdg(k_0) has positive imaginary part
%          %dkdg(k_{-1}) has negative imaginary part
%          [Y,I]=sort( -imag(dkdg));
%          inf_rts(1:2) = inf_rts(I);
%        else
%          %dkdg(k_0) has negative real part
%          %dkdg(k_{-1}) has positive real part
%          [Y,I]=sort( real(dkdg));
%          inf_rts(1:2) = inf_rts(I);
%        end
%      end
%    end
  
  
  
  %Use the Naewton-Raphson methods to find the first root for finite depth
  for nt = 1:1;
    x = inf_rts(nt);
    reps = 0;
    for ct = 1:100
      reps = reps + 1;
      %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      x_4 = C5*x^4;
      xd = x*fluid_depth(ht);
      tanch = tanh(xd);
      f = C0 + x*tanch*( C1 + x_4);
      df = tanch*( C1 + 5*x_4) + xd*sheck(xd)^2*( C1 + x_4);%sheck gives correct values for xd large and complex
      %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      if abs(f)<tol, break, end
      x = x - f / df;
    end
    if reps==100, 
      k(1,nt) = 0; 
    else
      k(1,nt) = x;
    end
  end
%    plot(real(k(1,1)), imag(k(1,1)), 'r.',real(k(1,2)), imag(k(1,2)), 'go',real(k(1,3)), imag(k(1,3)), 'bs')
%    pause
  
  %Zero viscosity evanescent modes
  k_evan_novisc = evanescent_modes(C0,C1r,C5,fluid_depth(ht),Dim);
  %Modify these using the tiny_steps approach
  if visc_rp > 0
    for nt = 1:Dim;
      x = k_evan_novisc(nt);
      for tt = 1:tinysteps_1 
        C1tt = C1r - C1i*tt;
        reps = 0;
        for ct = 1:100  %Newton's method
          reps = reps + 1;
          %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          x_4 = C5*x^4;
          xd = x*fluid_depth(ht);
          tanch = tanh(xd);
          f = C0 + x*tanch*( C1tt + x_4);
          df = tanch*( C1tt + 5*x_4) + xd*(1-tanch^2)*( C1tt + x_4);
          %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          if abs(f)<tol2, break, end
          x = x - f / df;
        end%ct
      end%tt
      k_evan(1,nt) = x;
      if reps==100, k_evan(1,nt) = 0; end
    end%nt
  else%viscosity is zero
    k_evan = k_evan_novisc;
  end
  
  if isempty(k_evan_novisc) == 0, k(1,2:Dim+1) = k_evan; end
  kn(ht,:) = k;
  
  %Checks for the correctitude of the roots.  
  %These error checks make the function take twice as long to complete, but they're worth it.
  %
  %Check the dispersion relation really is satisfied
  num_disp  = check_disp0(k, C0, C1, C5, fluid_depth(ht), tol,tol2);
  %Chesk that distances between points are not too small
  num_seps  = check_seps(k.');
  %Check if there are any NaN's in the roots
  num_nan   = sum(isnan(k));
  %Check if the roots are in teh correct quadrant, and roughly the correct position.
  num_pos   = check_pos0(k);
  %Check if there are any zeroes in the roots
  num_zero  = sum(k==0);
  %Find which roots satisfy the principle of the argument check;
  %This error checking routine makes the function take 50 times longer to complete
  if check_errs_flag ==1, 
    num_contour = check_principle_arg(k, C0, C1, C5, fluid_depth(ht));
  else
    num_contour = 0;
  end
  errs(ht,:) = [ num_nan  num_zero  num_disp  num_seps num_pos num_contour];

  clear k
end %for ht


return

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function k = evanescent_modes(C0,C1,C5,H,Dim)

% Find the evanescent modes of the Euler-Bernoulli (no viscosity) dispersion relation.
%
% k = evanescent_modes(C0,C1,C5,H)
%
% The relation is 
% C0 + k*tanh(k*H)*( C1 + C5*k^4) = 0
%
% Here C1 is real - no viscosity is included here.
%
% Rather than search for values of k, we search for values of g, where i*g=k*H; 
% The dispersion relation becomes  C0 - g*tan(g)*( C1/H + C5/H^5*g^4) = 0 
% but it seems better to solve
% C0/(g*tan(g)) - ( C1/H + C5/H^5*g^4) = 0.  First I use this form for all Dim modes and then, 
% for those roots that do not converge, I use
% C0*cos(g) - g*sin(g)*( C1/H + C5/H^5*g^4) = 0.
% then try
% C0*cos(g)/g - sin(g)*( C1/H + C5/H^5*g^4) = 0.
%
% If that doesn't work then an error message is written.
%
% Note that there is a zero for g \in [(n-1)*pi, n*pi], n a non-negative integer

tol = 1e-7;    %a critical parameter: avoid changing this: what works for one input, doesn't for another.
C1 = C1 / H;
C5 = C5 / H^5;
g(1:Dim) = complex(0);

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Try the first form of the dispersion relation
%
unconverged = find(g==0);
for nt = unconverged  %This run should converge for most of the roots
  left = (nt-1)*pi;   
  right = nt*pi;
  x = 0.5*left + 0.5*right;
  reps = 0;
  for ct = 1:100;  %Newton's method
    reps = reps+1;
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    tanx = tan(x);
    xtanx = x*tanx;
    f = C0/(xtanx) - (C1+C5*x^4);                      %The dispersion relation ...
    df = - C0*(tanx+x+x*tanx^2 )/(xtanx)^2 - 4*C5*x^3; %... and its derivative:
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    if abs(f)<tol , break, end
    x = x - f / df;
  end
  g(nt) = x;
  %Ascertain if the root can called converged
  if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0; end
end

%Identify unconverged roots
unconverged = find(g==0);

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Try the second form of the dispersion relation
%
for nt = unconverged  %for each of teh roots that have not converged
  left = (nt-1)*pi;
  right = nt*pi;
  x = right;          %Starting value
  reps = 0;
  for ct = 1:100;     %Newton's method
    reps = reps+1;
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    cx = cos(x);
    sx = sin(x);
    C5x4 = C5*x^4;
    T1 = (C1+C5x4);
    f = C0*cx - x*sx*T1;                   %The dispersion relation ...
    df = -sx*(C0 + C1 + 5*C5x4) - x*cx*T1; %... and its derivative:
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    if abs(f)<tol, break, end
    x = x - f / df;
  end
  g(nt) = x;
  %Ascertain if the root can called converged
  if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0;end
end
%Identify unconverged roots
unconverged = find(g==0);

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Try the third form of the dispersion relation
%
for nt = unconverged  %for each of the roots that have not converged
  left = (nt-1)*pi;
  right = nt*pi;
  x = right;          %Starting value
  reps = 0;
  for ct = 1:100;     %Newton's method
    reps = reps+1;
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    cx = cos(x);
    sx = sin(x);
    C5x4 = C5*x^4;
    T1 = (C1+C5x4);
    f = C0*cx/x - sx*T1;                   %The dispersion relation ...
    df = -sx*(C0 + 4*C5x4)/x - cx*( C0/x^2 + T1  ); %... and its derivative:
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    if abs(f)<tol, break, end
    x = x - f / df;
  end
  g(nt) = x;
  %Ascertain if the root can called converged
  if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0;end
end
%Identify unconverged roots
unconverged = find(g==0);


% Try the fourth form of the dispersion relation
%
for nt = unconverged  %for each of the roots that have not converged
  left = (nt-1)*pi;
  right = nt*pi;
  x = right;          %Starting value
  reps = 0;
  for ct = 1:100;     %Newton's method
    reps = reps+1;
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    tx = tan(x);
    C5x4 = C5*x^4;
    T1 = (C1+C5x4);
%      f = C0 - x*tx*T1;                   %The dispersion relation ...
%      df = -tx*(C1 + 5*C5x4)/x - x*T1*sec(x)^2; %... and its derivative:
    f = C0/T1 - x*tx;                   %The dispersion relation ...
    df = - C0/T1^2*4*C5x4/x -tx -x*sec(x)^2; %... and its derivative:
    %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    if abs(f)<tol, break, end
    x = x - f / df;
  end
  g(nt) = x;
  %Ascertain if the root can called converged
  if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0;end
end
%Identify unconverged roots
%unconverged = find(g==0);


if any(g==0), disp('there are unconverged roots in new_root_finder.m'); end
k = 1i*g/H;

return

function y = sheck(z)

%The value 710 is a bit of a bodge, but it seems to work. Much large and csch produces NaN

y = sech(  z.*(abs(real(z))<710 )   ).*(abs(real(z))<710 );

return

function y = check_disp0(k, C0, C1, C5, depth, tol,tol2)

N = length(k);

%  Find the error for the primary and complex roots
T1 = sum(abs( C0 + k(1:1) .*tanh(k(1:1) *depth).*( C1 + C5*k(1:1) .^4) ) > tol );

%Different forms of the dispersion relation are used to find different roots.
%Here the best disperion relation is found.
g = k(2:N)*depth / 1i;
g4 = g.^4;
cosg = cos(g);
sing = sin(g);
tang = sing./cosg;
gtang = g.* tang;
C1d = C1 / depth;
C5d5 = C5 / depth^5;
Cbrac = C1d + C5d5* g4;

%  5th form of the dispersion relation
S5 = abs( C0./Cbrac - gtang);       %This one seems to be a rather good estimator for the error.
T2 = sum(S5 > tol2);
%If necessary, use the other forms
if T2>0
  %  1st form of the dispersion relation
  S1 = abs( C0 - gtang.*Cbrac );
  %  2nd form of the dispersion relation
  S2 = abs( C0./gtang - Cbrac );
  %  3rd form of the dispersion relation
  S3 = abs( C0*cosg - g.*sing.*Cbrac);
  %  4th form of the dispersion relation
  S4 = abs( C0*cosg./g - sing.*Cbrac);
  T2 = sum(min([S1;S2;S3;S4;S5]) > tol2);
end

y = T1 + T2;

return

function num_seps = check_seps(kn)

N = length(kn);

dog = repmat( kn, [1 N]);
dog = abs(dog - dog.');
big23 = 99*max(max(dog));
big23 = big23 * eye(size(dog));
dog = dog + big23;

num_seps=sum(sum(dog<1e-8))/2;

return

function errs = check_pos0(kn)

N = length(kn);
Rmax = max(real(kn(2:N)));

%Check that the roots are in the correct quadrants
errs = sum(not( ( angle(kn)<=[ 0 zeros([1,N-1])]+pi/2 ) &  ...
                ( angle(kn)>=[ 0 zeros([1,N-1])]) ) );

%Check that they're in roughly the correct lpace: roots 1 and 2 to teh right of the evanescent roots.
if real(kn(1)) <= Rmax, errs = errs+1; end
%if real(kn(2)) <= Rmax, errs = errs+1; end

return

function [kn errs] = roots_rp(parameter_vector, visc_rp, Dim, h_vec, fluid_depth, check_errs_flag)

%Returns the dimensional roots of the Robinson-Palmer dispersion relation.
%Any values for period, visc_rp and h_vec may be used.

tol2 = 1e-5; %used in finding the shifted evanescent roots:critical, so don't change it 
tol = 1e-9; 
tinysteps_1 = 15;
%  tinysteps_2 = 15;
E=parameter_vector(4);%Pa
g=parameter_vector(5);%m/s^2
rho_w=parameter_vector(1);%kg/m^3
rho_ice=parameter_vector(2);%kg/m^3
nu=parameter_vector(3);
%  depth = 1000; %Used in the finite depth case
%  varpic = 5*(1/4)^.8 *real( (-1).^-.2 );


lh = length(h_vec);
kn=complex(zeros([lh, Dim+3]));
errs(1:lh,1:6)=complex(zeros([lh,6]));
  
    
for ht = 1:lh %For each thickness in the vector h_vec
  k = complex(zeros([1,Dim+3]));
  k_evan = k(1:Dim);
  om = sqrt(g*parameter_vector(6)); %2*pi/period;
  D = E* h_vec(ht)^3 / (12*(1-nu^2));
  C0 = -rho_w * om^2;
  C5 = D ;
%    L = (D/(rho_w*om^2))^0.2;
  C1r = rho_w * g  - om^2*rho_ice*h_vec(ht);
  C1i = 1i*visc_rp*om/tinysteps_1;
  C1 = C1r - 1i*visc_rp*om;
%    varpi = g / ( om^2 * L) - rho_ice * h_vec(ht) / ( rho_w*L); %Tim's non-dim param, but its useful here

  %starting values from the infinite depth case, with viscosity
  inf_rts = roots([C5 0 0 0 C1 C0]);
  [Y,I] = sort(2*pi*(angle(inf_rts)/pi <0 ) + angle(inf_rts));
  inf_rts = inf_rts(I);    %the roots, sorted by angle
  inf_rts(3) = -inf_rts(5);
  
  %For thickness in [0.001m 35m], period in [5s 35s] and for gamma in [0,5e3] the 
  %primary root has a smaller imaginary part than the shifted, 1st quadrant root.
  if imag(inf_rts(1)) > imag(inf_rts(2)) %Then they two roots need to be swapped
    inf_rts(1:2) = inf_rts(2:1);
  end

%    if abs(varpi-varpic)<0.01
%      %The dkdg method for distinguishing roo0ts doesn't always work.
%      %For varpi close to varpic use a tinysteps approach with a Newton-Raphson
%      
%      C1i = 1i*visc_rp*om/tinysteps_2;
%      
%      %starting values from the infinite depth case, zero viscosity
%      inf_rts = roots([C5 0 0 0 C1r C0]);
%      [Y,I] = sort(2*pi*(angle(inf_rts)/pi <0 ) + angle(inf_rts));
%      inf_rts = inf_rts(I);    %the roots, sorted by angle
%      inf_rts(3) = -inf_rts(5);
%  
%      for nt = 1:3;
%        for tt = 1:tinysteps_2
%          x = inf_rts(nt);
%          for ct = 1:100
%            CC1=C1r+C1i*tt;
%            %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%            f = C0 + CC1*x + C5*x^5;
%            df = CC1 + 5*C5*x^4;
%            %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%            if abs(f)<tol, break, end
%            x = x - f / df;
%          end
%          inf_rts(nt) = x;
%        end
%      end
%    else
%      %starting values from the infinite depth case, with viscosity
%      inf_rts = roots([C5 0 0 0 C1 C0]);
%      [Y,I] = sort(2*pi*(angle(inf_rts)/pi <0 ) + angle(inf_rts));
%      inf_rts = inf_rts(I);    %the roots, sorted by angle
%      inf_rts(3) = -inf_rts(5);
%     
%      
%      %It seems the derivative of the roots with respect viscosity provides a means for checking 
%      %which of the first quadrant roots is the primary and which is the complex, provided visc
%      dkdg = 1i*inf_rts(1:2).^2*L^2 ./ (4*L^5*inf_rts(1:2).^5 + 1);
%      %sort the values in inf_rts.
%      %for low viscosities angle alone distinguishes the two roots (and the dk/dg test doesn't work). 
%      if visc_rp*om > 0.1  %For larger viscosities use the dk/dg test
%        if varpi < varpic
%          %dkdg(k_0) has positive imaginary part
%          %dkdg(k_{-1}) has negative imaginary part
%          [Y,I]=sort( -imag(dkdg));
%          inf_rts(1:2) = inf_rts(I);
%        else
%          %dkdg(k_0) has negative real part
%          %dkdg(k_{-1}) has positive real part
%          [Y,I]=sort( real(dkdg));
%          inf_rts(1:2) = inf_rts(I);
%        end
%      end
%    end
  
  
  
  %Use the Naewton-Raphson methods to find the first three roots for finite depth
  for nt = 1:3;
    x = inf_rts(nt);
    reps = 0;
    for ct = 1:100
      reps = reps + 1;
      %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      x_4 = C5*x^4;
      xd = x*fluid_depth(ht);
      tanch = tanh(xd);
      f = C0 + x*tanch*( C1 + x_4);
      df = tanch*( C1 + 5*x_4) + xd*sheck(xd)^2*( C1 + x_4);%sheck gives correct values for xd large and complex
      %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
      if abs(f)<tol, break, end
      x = x - f / df;
    end
    if reps==100, 
      k(1,nt) = 0; 
    else
      k(1,nt) = x;
    end
  end
%    plot(real(k(1,1)), imag(k(1,1)), 'r.',real(k(1,2)), imag(k(1,2)), 'go',real(k(1,3)), imag(k(1,3)), 'bs')
%    pause
  
  %Zero viscosity evanescent modes
  k_evan_novisc = evanescent_modes(C0,C1r,C5,fluid_depth(ht),Dim);
  %Modify these using the tiny_steps approach
  if visc_rp > 0
    for nt = 1:Dim;
      x = k_evan_novisc(nt);
      for tt = 1:tinysteps_1 
        C1tt = C1r - C1i*tt;
        reps = 0;
        for ct = 1:100  %Newton's method
          reps = reps + 1;
          %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          x_4 = C5*x^4;
          xd = x*fluid_depth(ht);
          tanch = tanh(xd);
          f = C0 + x*tanch*( C1tt + x_4);
          df = tanch*( C1tt + 5*x_4) + xd*(1-tanch^2)*( C1tt + x_4);
          %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          if abs(f)<tol2, break, end
          x = x - f / df;
        end%ct
      end%tt
      k_evan(1,nt) = x;
      if reps==100, k_evan(1,nt) = 0; end
    end%nt
  else%viscosity is zero
    k_evan = k_evan_novisc;
  end
  
  if isempty(k_evan_novisc) == 0, k(1,4:Dim+3) = k_evan; end
  kn(ht,:) = k;
  
  %Checks for the correctitude of the roots.  
  %These error checks make the function take twice as long to complete, but they're worth it.
  %
  %Check the dispersion relation really is satisfied
  num_disp  = check_disp(k, C0, C1, C5, fluid_depth(ht), tol,tol2);
  %Chesk that distances between points are not too small
  num_seps  = check_seps(k.');
  %Check if there are any NaN's in the roots
  num_nan   = sum(isnan(k));
  %Check if the roots are in teh correct quadrant, and roughly the correct position.
  num_pos   = check_pos(k);
  %Check if there are any zeroes in the roots
  num_zero  = sum(k==0);
  %Find which roots satisfy the principle of the argument check;
  %This error checking routine makes the function take 50 times longer to complete
  if check_errs_flag ==1, 
    num_contour = check_principle_arg(k, C0, C1, C5, fluid_depth(ht));
  else
    num_contour = 0;
  end
  errs(ht,:) = [ num_nan  num_zero  num_disp  num_seps num_pos num_contour];

  clear k
end %for ht


return

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% function k = evanescent_modes(C0,C1,C5,H,Dim)
% 
% % Find the evanescent modes of the Euler-Bernoulli (no viscosity) dispersion relation.
% %
% % k = evanescent_modes(C0,C1,C5,H)
% %
% % The relation is 
% % C0 + k*tanh(k*H)*( C1 + C5*k^4) = 0
% %
% % Here C1 is real - no viscosity is included here.
% %
% % Rather than search for values of k, we search for values of g, where i*g=k*H; 
% % The dispersion relation becomes  C0 - g*tan(g)*( C1/H + C5/H^5*g^4) = 0 
% % but it seems better to solve
% % C0/(g*tan(g)) - ( C1/H + C5/H^5*g^4) = 0.  First I use this form for all Dim modes and then, 
% % for those roots that do not converge, I use
% % C0*cos(g) - g*sin(g)*( C1/H + C5/H^5*g^4) = 0.
% % then try
% % C0*cos(g)/g - sin(g)*( C1/H + C5/H^5*g^4) = 0.
% %
% % If that doesn't work then an error message is written.
% %
% % Note that there is a zero for g \in [(n-1)*pi, n*pi], n a non-negative integer
% 
% tol = 1e-7;    %a critical parameter: avoid changing this: what works for one input, doesn't for another.
% C1 = C1 / H;
% C5 = C5 / H^5;
% g(1:Dim) = complex(0);
% 
% % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% % Try the first form of the dispersion relation
% %
% unconverged = find(g==0);
% for nt = unconverged  %This run should converge for most of the roots
%   left = (nt-1)*pi;   
%   right = nt*pi;
%   x = 0.5*left + 0.5*right;
%   reps = 0;
%   for ct = 1:100;  %Newton's method
%     reps = reps+1;
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     tanx = tan(x);
%     xtanx = x*tanx;
%     f = C0/(xtanx) - (C1+C5*x^4);                      %The dispersion relation ...
%     df = - C0*(tanx+x+x*tanx^2 )/(xtanx)^2 - 4*C5*x^3; %... and its derivative:
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     if abs(f)<tol , break, end
%     x = x - f / df;
%   end
%   g(nt) = x;
%   %Ascertain if the root can called converged
%   if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0; end
% end
% 
% %Identify unconverged roots
% unconverged = find(g==0);
% 
% % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% % Try the second form of the dispersion relation
% %
% for nt = unconverged  %for each of teh roots that have not converged
%   left = (nt-1)*pi;
%   right = nt*pi;
%   x = right;          %Starting value
%   reps = 0;
%   for ct = 1:100;     %Newton's method
%     reps = reps+1;
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     cx = cos(x);
%     sx = sin(x);
%     C5x4 = C5*x^4;
%     T1 = (C1+C5x4);
%     f = C0*cx - x*sx*T1;                   %The dispersion relation ...
%     df = -sx*(C0 + C1 + 5*C5x4) - x*cx*T1; %... and its derivative:
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     if abs(f)<tol, break, end
%     x = x - f / df;
%   end
%   g(nt) = x;
%   %Ascertain if the root can called converged
%   if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0;end
% end
% %Identify unconverged roots
% unconverged = find(g==0);
% 
% % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% % Try the third form of the dispersion relation
% %
% for nt = unconverged  %for each of the roots that have not converged
%   left = (nt-1)*pi;
%   right = nt*pi;
%   x = right;          %Starting value
%   reps = 0;
%   for ct = 1:100;     %Newton's method
%     reps = reps+1;
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     cx = cos(x);
%     sx = sin(x);
%     C5x4 = C5*x^4;
%     T1 = (C1+C5x4);
%     f = C0*cx/x - sx*T1;                   %The dispersion relation ...
%     df = -sx*(C0 + 4*C5x4)/x - cx*( C0/x^2 + T1  ); %... and its derivative:
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     if abs(f)<tol, break, end
%     x = x - f / df;
%   end
%   g(nt) = x;
%   %Ascertain if the root can called converged
%   if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0;end
% end
% %Identify unconverged roots
% unconverged = find(g==0);
% 
% 
% % Try the fourth form of the dispersion relation
% %
% for nt = unconverged  %for each of the roots that have not converged
%   left = (nt-1)*pi;
%   right = nt*pi;
%   x = right;          %Starting value
%   reps = 0;
%   for ct = 1:100;     %Newton's method
%     reps = reps+1;
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     tx = tan(x);
%     C5x4 = C5*x^4;
%     T1 = (C1+C5x4);
% %      f = C0 - x*tx*T1;                   %The dispersion relation ...
% %      df = -tx*(C1 + 5*C5x4)/x - x*T1*sec(x)^2; %... and its derivative:
%     f = C0/T1 - x*tx;                   %The dispersion relation ...
%     df = - C0/T1^2*4*C5x4/x -tx -x*sec(x)^2; %... and its derivative:
%     %FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
%     if abs(f)<tol, break, end
%     x = x - f / df;
%   end
%   g(nt) = x;
%   %Ascertain if the root can called converged
%   if reps==100 || x<left || x>right || isnan(g(nt)), g(nt) = 0;end
% end
% %Identify unconverged roots
% %unconverged = find(g==0);
% 
% 
% if any(g==0), disp('there are unconverged roots in new_root_finder.m'); end
% k = 1i*g/H;
% 
% return

function y = check_disp(k, C0, C1, C5, depth, tol,tol2)

N = length(k);

%  Find the error for the primary and complex roots
T1 = sum(abs( C0 + k(1:3) .*tanh(k(1:3) *depth).*( C1 + C5*k(1:3) .^4) ) > tol );

%Different forms of the dispersion relation are used to find different roots.
%Here the best disperion relation is found.
g = k(4:N)*depth / 1i;
g4 = g.^4;
cosg = cos(g);
sing = sin(g);
tang = sing./cosg;
gtang = g.* tang;
C1d = C1 / depth;
C5d5 = C5 / depth^5;
Cbrac = C1d + C5d5* g4;

%  5th form of the dispersion relation
S5 = abs( C0./Cbrac - gtang);       %This one seems to be a rather good estimator for the error.
T2 = sum(S5 > tol2);
%If necessary, use the other forms
if T2>0
  %  1st form of the dispersion relation
  S1 = abs( C0 - gtang.*Cbrac );
  %  2nd form of the dispersion relation
  S2 = abs( C0./gtang - Cbrac );
  %  3rd form of the dispersion relation
  S3 = abs( C0*cosg - g.*sing.*Cbrac);
  %  4th form of the dispersion relation
  S4 = abs( C0*cosg./g - sing.*Cbrac);
  T2 = sum(min([S1;S2;S3;S4;S5]) > tol2);
end

y = T1 + T2;

return

function errs = check_pos(kn)


N = length(kn);
Rmax = max(real(kn(4:N)));

%Check that the roots are in the correct quadrants
errs = sum(not( ( angle(kn)<=[ 0 0 pi/2 zeros([1,N-3])]+pi/2 ) &  ...
                ( angle(kn)>=[ 0 0 pi/2 zeros([1,N-3])]) ) );

%Check that they're in roughly the correct lpace: roots 1 and 2 to teh right of the evanescent roots.
if real(kn(1)) <= Rmax, errs = errs+1; end
if real(kn(2)) <= Rmax, errs = errs+1; end


return

function [mu_0, mu_1] = mu_new_PWC_visc_wtd(parameter_vector, dim, ...
    root_vec_dr, mat_A, visc_rp)

% - 18.10.10 -> rewritten to include standard weights
% - modified from mu_new_PWC_May09 in `Uber Long'

freq = sqrt(parameter_vector(6)*parameter_vector(5));

gam = freq*visc_rp/parameter_vector(5)/parameter_vector(1);

%----------------------------------------------------------------%

[c_4, c_2, c_0] = quartic_eqn_new_PWC_visc_wtd(parameter_vector, dim, ...
    root_vec_dr, mat_A);

dum_alpha = parameter_vector(6)*parameter_vector(11);
dum_beta = parameter_vector(12);

% - 22.12.04 - %

v_sq = ((1-dum_alpha-1i*gam) / (c_4*dum_beta))...
     ...
     + ( c_0 / c_4 )...
     ...
     - ( ( c_2 / (2*c_4)) ^ 2 );
 

% - either v_sq is +ve or -ve and real - %
 
if v_sq > 0
    
    % - want the complex root in the upp left quad - %
    
    dummy = 1;
    
    v = sqrt(v_sq);
    
else
    
    % - want the upper half p.i. root or +ve real root - %
    
    dummy = 0;
    
    v = sqrt(-v_sq)*1i;
    
end
    
u = c_2 / (2*c_4);

lam_0 = -u + 1i*v;


% - either lam_0 is on the real line or it is in the upper half plane - %

mod_lam0 = abs(lam_0);
arg_lam0 = angle(lam_0); % \geq 0 \leq pi

lam_1 = -u - 1i*v;

% - either lam_1 is on the real line or it is in the lower half plane - %

mod_lam1 = abs(lam_1);
arg_lam1 = angle(lam_1); % \geq pi \leq 2*pi

mu_0 = sqrt(mod_lam0)*exp(1i*(arg_lam0/2));

mu_1 = sqrt(mod_lam1)*exp(1i*(arg_lam1/2) + 1i*dummy*pi);

return

function [c_4, c_2, c_0] = quartic_eqn_new_PWC_visc_wtd(parameter_vector, ...
    dim, root_vec_dr, A)

% - 18.10.10 -> rewritten to include standard weights

%----------------------------------------------------------------%

dum_H = parameter_vector(10);

%------------------------------------------------------------%

% one_vec = ones(dim,1);
% 
% W = diag(wt_vec_dr);
% 
% R = diag(root_vec_dr);
% 
% C = diag(cosh(root_vec_dr.*dum_H));
% 
% S = diag(sinh(root_vec_dr.*dum_H));
% 
% WSAWC = W*S*(A\W)*C; clear A W S C 
% 
% c_4 = 1;
% 
% c_2 = one_vec'*R*WSAWC*one_vec;
% 
% c_0 = one_vec'*(R^3)*WSAWC*one_vec;

one_vec = ones(dim,1);

R = diag(root_vec_dr);

T = diag(tanh(root_vec_dr.*dum_H));

c_4 = 1;

c_2 = one_vec'*R*T*inv(A)*one_vec;

c_0 = one_vec'*(R^3)*T*inv(A)*one_vec;

return

function C = mat_C_dr_PWC_wtd(parameter_vector, dim, root_vec, mu_0, mu_1, mat_A)

kappa = parameter_vector(6);            
dum_H = parameter_vector(10);
dum_beta = parameter_vector(12);

%--------------------------------------%
%This code replaced by Gareth 21/10/08
%  v_0 = vec_v_new_PWC_v2(parameter_vector, mu_0, dim, root_vec, section);
%  v_1 = vec_v_new_PWC_v2(parameter_vector, mu_1, dim, root_vec, section);

rv2 = root_vec.^2;
rth = root_vec.* tanh(root_vec*dum_H);
v_0 = -dum_beta*(mat_A\((rv2 + mu_0^2).*rth)) ;
v_1 = -dum_beta*(mat_A\((rv2 + mu_1^2).*rth)) ;

%--------------------------------------%


C = zeros(dim + 2);
C(1:dim,1:dim) = eye(dim);

%--------------------------------------%This code replace by Gareth 21/10/08
%  for loop = 1:dim
%      C(dim+1, loop) = root_vec(loop)*tanh(root_vec(loop)*dum_H)/kappa;
%      C(dim+2, loop) = -dum_beta * (root_vec(loop)^2) * C(dim+1, loop);
%  end
C(dim+1,1:dim) = rth/kappa;
C(dim+2,1:dim) = -dum_beta * rv2 .* C(dim+1, 1:dim).';
%--------------------------------------%

C(1:dim, dim+1) = v_0;
C(1:dim, dim+2) = v_1;
C(dim+1 , dim+1:dim+2) = ones(1,2);
C(dim+2, dim+1) = -dum_beta*(mu_0^2);
C(dim+2, dim+2) = -dum_beta*(mu_1^2);

return

% WATER ICE TRANSITION

function [Rm, Tm, Rp, Tp]  = ...
    CalcWaterIceTransition_wtd(parameter_vector, Dim, DimG, DimI0, DimI,...
    k0, kk, mu_vec, C_mat, mat_A0, mat_A, Geom, yn_cmpx)

% - 18.10.10 -> rewritten to include standard weights

DD=Geom(1); dr=Geom(2); hh=Geom(3); al=Geom(4); be=Geom(5);

theta=0;

% [Rm, Tm, Rp, Tp]  = ...
%     CalcWaterIceTransition_v3(parameter_vector, Dim, DimG, ...
%     k0, kk, wt_0, wt, mu_vec, DD, dr, hh, al, be, PM)

% - LB 27.07.10 - updated from CalcWaterIceTransition_v3.m for WIFAR
%               - for use in multiple scattering probs      

Tol_vec(1) = 1e-16; % - Real root error - %
Tol_vec(2) = 1e-16; % - Imag root error (init) - %
%Tol_vec(3) = 1e-4; % - Tol on const soln - %

% - has to be done for each incident angle, the value u_n is an integer
% that defines which of the values of angles we take

% al = (parameter_vector(2)/parameter_vector(1))*DD;
% be = parameter_vector(4)*(DD^3)/...
%     (12*parameter_vector(1)*(1-(parameter_vector(3)^2))*parameter_vector(5));

parameter_vector = [parameter_vector(1:6), hh, dr, DD, hh-dr, al, be];

clear al be

% - calc flat bed fs onto uni ice

% - S = Scat Matrix; P = Trans Matrix

% - u = periodicity in y-direction (retained)
% - PM defines whether the ice-cover stretches to 1=+infty -1=-infty
%PM=1;
% - DD is ice thickness; dr is draught; hh is bed

%kk = zeros(Dim,1); k0 = zeros(Dim,1); wt = zeros(Dim,1); wt_0 =
%zeros(Dim,1); 
%for loop_Dim = 1:Dim
%kk(loop_Dim) = GetRootsMMA_PWC(parameter_vector, loop_Dim, Tol_vec);
%wt(loop_Dim) = weight_PWC(parameter_vector, kk(loop_Dim));
%k0(loop_Dim) = GetRootsMMA_FS_PWC(parameter_vector, loop_Dim, Tol_vec); 
%wt_0(loop_Dim) = weight_0_PWC(parameter_vector, k0(loop_Dim));
%end

%[mu_vec(1), mu_vec(2)] = mu_new_PWC(parameter_vector, Dim, kk, wt);

u = k0(1)*sin(theta);

til_k0 = sqrt(k0.^2 - u^2); til_kk = sqrt(kk.^2 - u^2);
til_mu0 = sqrt(mu_vec(1)^2 - u^2); 
til_mu1 = sqrt(mu_vec(2)^2 - u^2); til_mu1 = -til_mu1;

% - Required Matrices - %

%mat_A = matrix_A_PWC(parameter_vector, Dim, kk, wt, parameter_vector(10));
%mat_A0 = matrix_A_PWC(parameter_vector, Dim, k0, wt_0, parameter_vector(7));

%% - version: u in vertical modes

mat_W = mat_A;

mat_W0 = fn_jump_W(parameter_vector, Dim, DimG, k0, kk);

mat_W0_T = conj(mat_W0'); mat_W_T = conj(mat_W');

% --------------------------------------------------------- %

% - Calc the amps of the complex waves in terms of the prop waves

%til_ps = 1-parameter_vector(3);

% - Nb. in is from +infty, out is towards +infty

Bend_kk_in = C_mat(Dim+2,1:Dim); % + parameter_vector(12)*til_ps*(u^2)*C_mat(Dim+1,1:Dim);
Bend_mu0_in = C_mat(Dim+2,Dim+1); % + parameter_vector(12)*til_ps*(u^2)*C_mat(Dim+1,Dim+1);
Bend_mu1_in = C_mat(Dim+2,Dim+2); % + parameter_vector(12)*til_ps*(u^2)*C_mat(Dim+1,Dim+2);

Shear_kk_in= -transpose(til_kk).*C_mat(Dim+2,1:Dim); % - parameter_vector(12)*til_ps*(u^2)*C_mat(Dim+1,1:Dim) );
Shear_mu0_in = -til_mu0*C_mat(Dim+2,Dim+1); % - parameter_vector(12)*til_ps*(u^2)*C_mat(Dim+1,Dim+1) );
Shear_mu1_in = -til_mu1*C_mat(Dim+2,Dim+2); % - parameter_vector(12)*til_ps*(u^2)*C_mat(Dim+1,Dim+2) );

Mat_muOut = [[Bend_mu0_in; -Shear_mu0_in] [Bend_mu1_in; -Shear_mu1_in]];
Mat_In = [[Bend_kk_in; Shear_kk_in] [Bend_mu0_in; Shear_mu0_in] [Bend_mu1_in; Shear_mu1_in]];
vec_kkOut = [Bend_kk_in; -Shear_kk_in];

% vec_muIn = -Mat_muOut\vec_kkIn; 
MatBS_Out = -Mat_muOut\vec_kkOut;
% vec_mu = -Mat_muOut\Mat_muIn;
MatBS_In = -Mat_muOut\Mat_In;
clear Mat_muOut vec_kkOut

% -------------------------------------- %

phi0In = eye(Dim); phi0Out = eye(Dim);
phi0In_dx = 1i*diag(til_k0); phi0Out_dx = -1i*diag(til_k0);

%psiIn = zeros(1:Dim,1:Dim+2); psiIn_dx = psiIn;

psiIn = C_mat(1:Dim,1:Dim+2) + C_mat(1:Dim,Dim+1:Dim+2)*MatBS_In; 
%psiIn(:,Dim+[1:2]) = C_mat(1:Dim,Dim+1:Dim+2)*(eye(2)+vec_mu);
psiOut = C_mat(1:Dim,1:Dim) + C_mat(1:Dim,Dim+1:Dim+2)*MatBS_Out;
psiIn_dx = -1i*C_mat(1:Dim,1:Dim+2)*diag([til_kk;til_mu0;til_mu1]) +...
    1i.*C_mat(1:Dim,Dim+1:Dim+2)*diag([til_mu0, til_mu1])*MatBS_In;
psiOut_dx = 1i*C_mat(1:Dim,1:Dim)*diag(til_kk) +...
    1i.*C_mat(1:Dim,Dim+1:Dim+2)*diag([til_mu0, til_mu1])*MatBS_Out;

%% - 03.07.09: non-square mats mat_W and mat_W0

% - following are Dim x DimG:

MatB_u0 = (phi0Out_dx\(mat_A0\mat_W0)); 
MatB_u = (psiOut_dx\(mat_A\mat_W)); 

% - following are Dim x Dim+2:

MatB_I = -psiOut_dx\psiIn_dx;

% - following are Dim x Dim:

MatB_I0 = -phi0Out_dx\phi0In_dx;

% - following are DimG x Dim:

MatA_I0 = mat_W0_T*(phi0In + phi0Out*MatB_I0);

% - following are DimG x Dim+2:

MatA_I = mat_W_T*(psiIn + psiOut*MatB_I);

% - following are DimG x DimG:

MatA_u0 = mat_W0_T*phi0Out*MatB_u0;
MatA_u = mat_W_T*psiOut*MatB_u;

% - find u in terms of the inc amps (DimG x Dim)

MatU_I0 = (MatA_u - MatA_u0)\MatA_I0;

% - find u in terms of the inc amps (DimG x Dim+2)

MatU_I = (MatA_u0 - MatA_u)\MatA_I;

% - Scat mat
v1 = 1:Dim; v2=Dim+1:Dim+2;

Rm = zeros(Dim,Dim); Rp = zeros(Dim+2,Dim+2); 
Tm = zeros(Dim+2,Dim); Tp = zeros(Dim,Dim+2);

Rm = MatB_I0 + MatB_u0*MatU_I0; 
Tp =           MatB_u0*MatU_I; 
Rp(v1,:) = MatB_I + MatB_u*MatU_I;
Tm(v1,:) =          MatB_u*MatU_I0;

Rp(v2,:) = MatBS_In + MatBS_Out*Rp(v1,:);
Tm(v2,:) =            MatBS_Out*Tm(v1,:);

% - For Output -> swap dimension

if DimI>Dim
 disp('error -> Output dimension greater than input')
elseif yn_cmpx == 1
 v1 = 1:DimI0; v2 = [1:DimI,Dim+1,Dim+2];
 
 Rm=Rm(v1,v1); Tp=Tp(v1,v2); Tm=Tm(v2,v1); Rp=Rp(v2,v2);   
elseif yn_cmpx == 0
 v1 = 1:DimI0; v2 = 1:DimI; 
 
 Rm=Rm(v1,v1); Tp=Tp(v1,v2); Tm=Tm(v2,v1); Rp=Rp(v2,v2);   
end
      
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUBFUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function WW = fn_jump_W(pv, dimN, dimM, rt_0, rt)
    
% - V is NxM - %

water_depth_j = pv(10); water_diff_i = pv(8);

WW = zeros(dimN,dimM); 

for loopN = 1:dimN
    
    for loopM=1:dimM
        
        WW(loopN,loopM) = ...
            ...
            fn_wip(rt_0(loopN,:), rt(loopM,:), ...
            ...
            water_depth_j, water_diff_i);
        
    end
    
end

return

function cosh_cosh = fn_wip(root_i, root_j, water_depth_j, water_diff_i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Luke
%Abbreviated by Gareth
% - Rewritten by Luke to inc the weights sech(k_i H)
% - 3rd Sept 08 - Rewritten by Luke to deal with subtractive cancelation
% - Rewritten by Luke to avoid redundant calls
% - Rewritten by Gareth to deal with v.short waves and a v.small difference in roots

%Evaluates integrals of the form
%int_{-h}^{-d} cosh( root_i(it)*(z+h) ) cosh( root_j(jt)*(z+h) ) dz
%where z=-h is the seafloor, and -d is the underside of the ice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 1e-8; % SET THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tol3 = -100; % SET THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

it = 1; jt = 1;

if abs(root_i(it,1) - root_j(jt,1))<tol %% when the 2 wavenumbers get close
    eps = root_i(it,1)-root_j(jt,1);
    arg = 2*water_depth_j*root_j(jt,1);
    if abs(arg) >708
        cosh_cosh(it,jt) = 1 / (root_i(it,1) + root_j(jt,1));
    else
        term = ( arg + sinh(arg) )/(4*root_j(jt,1)); clear arg
        cosh_cosh(it,jt) = term;
        loop_eps = 0;
        %       tol2 = 1e-15; %% SET THIS!!!!!!!!!!
        while loop_eps<2 %abs((eps^loop_eps)*term)<tol2
            loop_eps=loop_eps+1;
            term = jump_exp(root_j(jt,1), water_depth_j, loop_eps);
            cosh_cosh(it,jt) = cosh_cosh(it,jt) + (eps^loop_eps)*term;
        end
        cosh_cosh(it,jt) = sech(root_i(it,1)*(water_depth_j+water_diff_i))*...
            sech(root_j(jt,1)*water_depth_j)*cosh_cosh(it,jt);
    end
else %% if roots are distinct & taking into account the weighting
    water_depth_i = water_depth_j + water_diff_i;
    %for root_i large and deep water tanh(root_i(it,1)*water_depth_i) ~ 1: use the else stuff.
    %The condition derives from the leading error term in the else stuff approximation.
    if real(root_i(it,1)*( water_diff_i-2*water_depth_i)) >  tol3
        sinh_i = tanh(root_i(it,1)*water_depth_i)*cosh(root_i(it,1)*water_diff_i) - sinh(root_i(it,1)*water_diff_i);
        cosh_i = cosh(root_i(it,1)*water_diff_i) - tanh(root_i(it,1)*water_depth_i)*sinh(root_i(it,1)*water_diff_i);
        tanh_j = tanh(root_j(jt,1)*water_depth_j);
        quot = root_i(it,1)^2 - root_j(jt,1)^2;
        quot = 1 / quot;
        cosh_cosh(it,jt) = quot*(root_i(it,1)*sinh_i - root_j(jt,1)*tanh_j*cosh_i);
    else
        sinh_i = exp(- root_i(it,1)*water_diff_i);
        cosh_i = sinh_i;
        tanh_j = tanh(root_j(jt,1)*water_depth_j);
        quot = root_i(it,1)^2 - root_j(jt,1)^2;
        quot = 1 / quot;
        cosh_cosh(it,jt) = quot*(root_i(it,1)*sinh_i - root_j(jt,1)*tanh_j*cosh_i);
    end
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunction to calculate the terms in the Taylor exp of the above when
%% the roots are close

function T = jump_exp(k, H, N)

% display(['in subfn', N])

arg = 2*k*H;

if N==1
    T = ( arg*cosh(arg) - sinh(arg) )/(8*(k^2));
elseif N==2
    T = ( (arg^3/2) - 3*arg*cosh(arg) + 3*(1+(arg^2)/2)*sinh(arg) )/(48*(k^3));
elseif N==3
    T = ( arg*(3+(arg^2)/2)*cosh(arg) - 3*(1+(arg^2)/2)*sinh(arg) )/(96*(k^4));
elseif N==4
    T = ( (arg^5)/(2^4) - (5/2)*arg*(3+(arg^2)/2)*cosh(arg) + (5/2)*(3+(3/2)*(arg^2)+(1/8)*(arg^4))*sinh(arg) )/(480*(k^5));
elseif N==5
    T = ( arg*(15+(5/2)*(arg^2)+(1/8)*(arg^4))*cosh(arg) - 5*(3+(3/2)*(arg^2)+(1/8)*(arg^4))*sinh(arg) )/(1920*(k^6));
else
    display('need more terms in taylor expansion')
end

return

% - NUMERICAL INTEGRATION - %

% approximate int_{a}^{b}f(x)dx = fac*wts*f(x_abs)

function [wts,x_abs,fac] = fn_NumInt(Nint,x0,x1)

%%% CHEBYSHEV %%%

% fac=(x1-x0)/2;
% 
% [x,w]=OP_numint_chebyshev(Nint);
% 
% x_abs = fac*(x+1) + x0;
% 
% wts = (w.*sqrt(1-x.^2)).';

%%% LEGENDRE %%%

fac=(x1-x0)/2;

[x,wts]=OP_numint_legendre(Nint,[-1,1]);

x_abs = fac*(x+1) + x0;

wts = wts.';

return

% this module provides the subroutine GEN_gauleg, which computes
% N-vectors x & w of abscissae and weights
% to use for Gauss-Legendre integration
% given N and the endpoints x1 & x2.
% CALL [x,w]=GEN_gauleg(x1,x2,N)function

function [x,w]=OP_numint_legendre(N,x1x2)

if nargin==1
  x1=-1;
  x2=1;
else
  x1=x1x2(1);
  x2=x1x2(2);
end

ndp=14;
EPS=5/10^(1+ndp);
MAXIT=10;

M=floor((N+1)/2); % The roots are symmetric in the interval,
          % so we only have to ﬁnd half of them.
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);

% Initial approximations to the roots:
zz=cos(pi*((1:M)'-0.25)/(N+0.5));
ppvec=0*zz;
zzvec=ppvec;
x=zeros(N,1);
w=x;

for jj=1:M
  unfinished=1;
  z=zz(jj);
  % Newton’s method carried out individually on the roots.
  for its=1:MAXIT
    p1=1;
    p2=0;
    % Loop up the recurrence relation to get
    % the Legendre polynomials evaluated at z
    if unfinished
      for j=1:N
         p3=p2;
         p2=p1;
         p1=((2*j-1)*z*p2-(j-1)*p3)/j;
      end
      % p1 now contains the desired Legendre polynomials.
      % We next compute pp, the derivatives, by a standard relation
      % involving also p2, the polynomials of one lower order.
      pp=N*(z*p1-p2)/(z*z-1);
      z1=z;
      z=z1-p1/pp;
      unfinished=(abs(z-z1) > EPS);
    else
      zzvec(jj)=z;
      ppvec(jj)=pp;
      break
    end
  end
  if (its == MAXIT+1)
    disp('too many iterations in GEN_gauleg')
  end
end

%OUTPUTS
x(1:M)=xm-xl*zzvec;        % Scale the root to the desired interval
x(N:-1:N-M+1)=xm+xl*zzvec; % Put in its symmetric counterpart.
w(1:M)=2*xl./((1-zzvec.^2).*ppvec.^2);  % Compute the weight
w(N:-1:N-M+1)=w(1:M);                   % and its symmetric counterpart

return

