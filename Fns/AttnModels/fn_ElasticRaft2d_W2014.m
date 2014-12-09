% function out = ...
%  fn_ElasticRaft2d_WP2009(fortyp,lam0,Param,outputs,SURGE,LONG,COMM,Ens_size,DO_PLOT,col)
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
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam_vec = a SCALAR of forcing fortyp

%% Param is structure with attributes and default values (Param_def_Oceanide.m):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Np: 1                         no of disks                   [-]
%%   thickness: 0.033000000000000         thickness                     [m]
%%           g: 9.810000000000000         gravity                       [m s^{-2}]
%%       rho_0: 1025                      water density                 [kg m^{-3}]
%%         rho: 5.590909090909091e+02     floe density                  [kg m^{-3}]
%%          nu: 0.300000000000000         Poisson's ratio               [-]
%%           E: 1.000000000000000e+10     Youngs modulus                [Pa]
%%       draft: 0.018000000000000         draft, rho/rho_0*thickness    [m]
%%           D: 3.290934065934066e+04     flexural rigidity             [Pa m^3]
%%        beta: 3.272851561059214         scaled rigidity, D/rho/g      [m^4]
%%         bed: 3.100000000000000         water depth                   [m]
%%  MIZ_length: 5                         length of MIZ                 [m]
%%   floe_diam: 0.990000000000000         floe diameter                 [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputs:
%'wavenumber'              open water wavenumber                  k0    
%'reflected energy'        R^2 for 1 floe                         RR    
%'transmitted energy'      T^2 for 1 floe                         TT    

%% HEAVE & PITCH:
%% eta=[heave RAO]+[pitch RAO]*x for -D/2<x<D/2 (ie D=floe_length)
%%   when amplitude=1m
%% [heave RAO] = \int_{-D/2}^{D/2}(eta)dx/I0       ~ (fac/I1)*(wts.*eta)
%%    I0       = \int_{-D/2}^{D/2}(1)dx=D
%% [pitch RAO] = \int_{-D/2}^{D/2}(eta*x)dx/I1     ~ (fac/I1)*(wts.*x.*eta)
%%    I1       = \int_{-D/2}^{D/2}(x^2)dx=D^3/12

%% 'heave'                 [heave RAO] (see above)
%% 'pitch/k'               P=[pitch RAO]/k0 => pitch=[pitch RAO]*inc_amp=P*(ka) ie P=pitch/steepness
%% 'pitch-ang'             180*P/pi

%% SURGE:
%% 'surge-nonorm'          surge RAO (m for 1m amplitude)      abs(u_sg)
%% 'surge-norm'            surge normalised by ??????????      abs(u_sg*tanh(k0*depth))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SURGE=0/1 surge on/off
% LONG=0/1  long floe approx on/off
% COMM = comments on or off
% Ens_size = how many members of ensemble
%            only necessary when long floe limit off


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
 fn_2dFloe_WP2009(fortyp,lam0,Param,outputs,SURGE,LONG,COMM,Ens_size,DO_PLOT,col)

%Model used
smod  = 'Williams & Porter (2009) model';

%% Inputs & prelims

if ~exist('DO_PLOT','var'); DO_PLOT=0; end
if ~exist('col','var');     col='r'; end

if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('lam0','var');   lam0=1/2.5; end

if ~exist('LONG','var'); LONG=0; end
%if LONG==0
%   LONG  = 1;
%   disp(['Only LONG=1 activated so far in ',smod,':'])
%   disp(['Setting LONG=1'])
%end

if ~exist('COMM','var'); COMM=1; end

if ~exist('SURGE','var'); SURGE=0; end
if SURGE==1
   SURGE = 0;
   disp(['Surge not yet possible in ',smod,':'])
   disp(['Setting SURGE=0'])
end


%if ~exist('RIGID','var'); RIGID=10; end

if ~exist('Param','var');
   Param = ParamDef_Oceanide(5); 
 end

if ~LONG; if ~exist('Ens_size','var'); Ens_size = 1; end; end

if ~exist('outputs','var'); 
 outputs='transmitted energy heave pitch/k surge-norm'; end

depth = Param.bed;

Forcing = Force_def(Param.g(1), depth, fortyp, lam0);


thickness = Param.thickness;

floe_length = Param.floe_diam;

%Vert_Dim = Param.Ndtm;
%Inc = eye(Vert_Dim,1);
%DimG = Vert_Dim; clear Dims

draught = Param.draft; al = draught;

be = Param.beta;

fq = (2*pi*Forcing.f)^2/Param.g;

parameter_vector = [Param.rho_0, Param.rho, Param.nu, Param.E, ...
    Param.g, fq]; 
   

%% %-------------------------------------------------------------------------%
%% % Surge
%% %
%% % Kinematic condition
%% %
%% % phi = kappa*u*cos(theta)
%% %
%% % Dynamic condition
%% %
%% % (S-kappa*M-i*A)*u =
%% %     rho_0*int_{-d}^{0}{phi(L/2,z)-phi(-L/2,z)}dz
%% 
%% modMass = floe_length*draught;
%% damp = 0;
%% spring = 0;
%% % adjust these for 2d problem is neq 0
%% modDamp = 0; %2*Forcing.f*damp/Param.g/Param.rho_0;
%% modSpring = 0; %spring/Param.rho_0/pi;
%% clear damp spring
%% 
%% %-------------------------------------------------------------------------%

%visc_rp = 0; %1e-1; %1e-1; % 

%% %% FREE-SURF ROOTS
%% 
%% %Roots0 = roots_rp0(parameter_vector, 0, Vert_Dim-1, 0, depth, 0);
%% %Roots0=Roots0.';
%% 
%% Roots0 = zeros(Vert_Dim,1); 
%% 
%% Roots0(1)   = fn_RealRoot([depth*parameter_vector(6)], ...
%%           'fn_ReDispRel_water', 'fn_UppLimReal_water', 1e-16)/depth; 
%% 
%% for loop_Dim = 2:Vert_Dim
%%  Roots0(loop_Dim) = 1i*fn_ImagRoot_water(loop_Dim-1, depth, ...
%%   parameter_vector(6), 1e-16); 
%% end
%% 
%% mat_A0 = wtd_cosh_cosh(Vert_Dim, Roots0, depth);    
%% if SURGE; mat_Q0 = fn_Qmat(depth,draught,Roots0); 
%% else mat_Q0=0*mat_A0(:,1); end

%% Begin calculations:

if COMM
  cprintf([0.3,0.3,0.3],'------------------------------------------\n')
  cprintf([0.3,0.3,0.3],'----------   START: 2d Floe    -----------\n')
end

if COMM
 cprintf([0.3,0.3,0.3],['>>> ' fortyp ' = ' num2str(lam0) '\n'])
 if ~LONG 
  if SURGE
   cprintf([0.3,0.3,0.3],['>>> with surge \n'])
  else
   cprintf([0.3,0.3,0.3],['>>> without surge \n'])
  end
  cprintf([0.3,0.3,0.3],['>>> floe length = ' num2str(floe_length) '\n'])
  cprintf([0.3,0.3,0.3],['>>> ensemble = ' num2str(Ens_size) '\n'])
 else
  cprintf([0.3,0.3,0.3],'>>> long floe limit\n')
 end
 cprintf([0.3,0.3,0.3],['>>> ' num2str(thickness) ' thick\n'])
 cprintf([0.3,0.3,0.3],['>>> rigidity = ' sprintf('%0.5g',Param.E) '\n'])
end
     
%clear Param
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ONE INTERFACE 

%[Rm0,Tm0,Rp0,Tp0,Sm0_sg,Sp0_sg,Roots,c_vec] = ...
%            fn_WaterIce(parameter_vector, Vert_Dim, DimG, Roots0, ... 
%            mat_A0, mat_Q0, thickness, depth, draught, al, be, SURGE);

%%inputs for SUB_RTstep_Gal_manymodes.m 
hh          = [0,Param.thickness];
period      = 1/Forcing.f; 
H1          = depth;
H2          = depth-draught;
phys_vars   = {period,0,H2};
bc          = 1;%%free edge conditions;
if LONG==1
   M_inner  = 1;
else
   M_inner  = 50;
end
MM = [1 M_inner];
NN = [50 1000];%% # Gegenbauer poly's, # vertical modes

INC_SUB  = 1;%%Include draught
EE       = [0,Param.E;
            0,Param.rho;
            0,Param.nu];
rho_wtr  = Param.rho_0;

[Rp,Tp,Rm,Tm,Smat,Y_Gal] = SUB_RTstep_Gal_manymodes(...
			 phys_vars,hh,bc,MM,NN,INC_SUB,EE,rho_wtr);
%%water wavenumbers;
Roots0   = Y_Gal{1}{1}; 
k0       = Roots0(1);%%real root
if COMM
   cprintf([0.3,0.3,0.3],['>>> ' int2str(NN(1)) ' Gegenbauer polynomials\n'])
   cprintf([0.3,0.3,0.3],['>>> ' int2str(NN(2)) ' vertical modes\n'])
   cprintf([0.3,0.3,0.3],['>>> water lam0/k0 = ' num2str(2*pi/Roots0(1)) '/' num2str(Roots0(1)) '\n'])
end

%%ice wavenumbers;
Roots = Y_Gal{1}{2}; 
k_ice = Roots(1);%%real root
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
if COMM; cprintf([0.3,0.3,0.3],['>>> ice lam/k   = ' ...
  num2str(2*pi/Roots(1)) '/' num2str(Roots(1)) '\n']); end           
        
if LONG %%% LONG FLOE LIMIT %%%
 
   r11 = Rp(1,1); 
   %alpha = -2*log(1-abs(r11)^2);
 
   TT = (1-abs(r11)^2)^2;
   RR = [];

else %%% NO LONG FLOE LIMIT %%% 
 
   sd = min(floe_length,pi/Roots(1));         
  
   for loop=1:Ens_size 
      
      if Ens_size==1
       dum_fl = floe_length;
      else
       dum_fl = floe_length + 2*sd*(rand - 0.5);
      end%% use transition matrices to get floe results

      TM0         = {Rp,Tp,Rm,Tm};
      [R,B0,B1,T] = FB_2steps_Tmatrices(TM0,[],Roots*floe_length);%%amp's in potential
      if Ens_size==1
         %%displacement amplitudes in floe (for unit amplitude from left);
         gam1        = Roots0(1:M_inner);
         gam2        = Roots(1:M_inner);
         disp_fac1   = gam1.*tanh(gam1*H1);
         disp_fac2   = gam2.*tanh(gam2*H2);
         b0_w        = B0.*disp_fac2/disp_fac1(1);%%waves from left edge in expansion for displacement
         b1_w        = B1.*disp_fac2/disp_fac1(1);%%waves from right edge in expansion for displacement
      end

      %[Rm,Tm,Rp,Tp,Sm_sg,Sp_sg] = ...
      %          fn_IndFloe(dum_fl, Rm0, Tm0, Rp0, Tp0, Sm0_sg, Sp0_sg, Roots);
       
      if SURGE
         %%not done
      else
         u_sg  = 0;
      end
      
      r_vec(loop) = abs(R(1))^2;
      t_vec(loop) = abs(T(1))^2;
   end
   
   Rm = mean(r_vec);
   Tm = mean(t_vec);
   TT = Tm;
   RR = Rm;
   
   if abs(RR+TT-1)>1e-3
    cprintf('magenta',['>>> warning R+T=' num2str(RR+TT) '\n'])
   end

end

%% PLOT:

if or(DO_PLOT,or(~isempty(strfind(outputs,'heave')),...
  ~isempty(strfind(outputs,'pitch'))))
 
   x_res          = 501;
   [wts,xx,fac]   = fn_NumInt(x_res,-floe_length/2,floe_length/2);

   ExL      = exp(1i*(xx+floe_length/2)*Roots(1:M_inner).');
   ExR      = flipud(ExL);
   eta      = ExL*b0_w+ExR*b1_w;
   eta_inc  = exp(1i*Roots0(1)*xx);
 
   %%% Plot
   if DO_PLOT
      figure(DO_PLOT)
      
      h1 = subplot(2,1,1); hold on; set(h1,'box','on')
      h2 = subplot(2,1,2); hold on; set(h2,'box','on')
      
      %%incident wave
      plot(h1,xx_ext,real(exp(1i*Roots0(1)*(xx_ext+floe_length/2))),'k:')
      %%floe profile
      plot(h1,xx,real(eta),col)
      ylabel(h1,'Re(\eta)','fontsize',14)
      title(h1,['floe profile (' col ') and incident wave (k:)'],'fontsize',14)

      %%incident wave
      plot(h2,xx_ext,imag(exp(1i*Roots0(1)*(xx_ext+floe_length/2))),'k:')
      %%floe profile
      plot(h2,xx,imag(eta),col)
      xlabel(h2,'x','fontsize',14); ylabel(h2,'Im(\eta)','fontsize',14)
   end
 
end

 
%% OUTPUTS & FINISH 

out_str = ' ''dummy'' '; out_val = ' 0 ';

if strfind(outputs,'wavenumber')
 out_str = [out_str '; ''k0'' '];
 out_val = [out_val '; Roots0(1)'];
end

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

if strfind(outputs,'pitch/k')
 out_str = [out_str '; ''pitch/k'' '];
 out_val = [out_val '; abs((12*fac/floe_length^3)*wts*(eta.*xx))/Roots0(1)'];
end

if strfind(outputs,'pitch-ang')
 out_str = [out_str '; ''pitch-ang'' '];
 out_val = [out_val ...
  '; 180*abs((12*fac/floe_length^3)*wts*(eta.*xx))/Roots0(1)/pi'];
end

if strfind(outputs,'surge-norm')
 out_str = [out_str '; ''surge'' '];
 out_val = [out_val '; abs(u_sg*tanh(Roots0(1)*depth))'];
end

if strfind(outputs,'surge-nonorm')
 out_str = [out_str '; ''surge'' '];
 out_val = [out_val '; abs(u_sg)'];
end

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});']) 
out(1)=[];

if COMM
 if strfind(outputs,'reflected energy')
  cprintf('blue',['>>>> reflected energy  : ' num2str(RR) '\n'])
 end
 if strfind(outputs,'transmitted energy')
  cprintf('blue',['>>>> transmitted energy: ' num2str(TT) '\n'])  
 end
 if strfind(outputs,'heave')
  cprintf('blue',['>>>> heave: ' ...
   num2str(abs((fac/floe_length)*wts*eta)) ...
   ' (inc wave: ' num2str(abs((fac/floe_length)*wts*eta_inc)) ') \n'])
 end
 if strfind(outputs,'pitch')
  cprintf('blue',['>>>> pitch/k: ' ...
   num2str(abs((12*fac/floe_length^3)*wts*(eta.*xx))/Roots0(1)) ...
   ' (inc wave: ' ...
   num2str(abs((12*fac/floe_length^3)*wts*(eta_inc.*xx))/Roots0(1)) ...
   ') \n'])
 end 
 if strfind(outputs,'surge-norm')
  cprintf('blue',['>>>> surge*tanh(kh): ' ...
   num2str(abs(u_sg.*tanh(Roots0(1)*depth))) ...
   '\n'])
   %' (eccentricity=' num2str(coth(Roots0(1)*depth)) ') \n'])
 end 
 if strfind(outputs,'surge-nonorm')
  cprintf('blue',['>>>> surge: ' ...
   num2str(abs(u_sg)) ...
   '\n'])
   %' (eccentricity=' num2str(coth(Roots0(1)*depth)) ') \n'])
 end 
 cprintf([0.3,0.3,0.3],'-----------    END: 2d Floe   ------------\n')
 cprintf([0.3,0.3,0.3],'------------------------------------------\n')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUBFUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
