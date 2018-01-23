function [alp,WW,tau_init,T_relax] = fn_ode_drag_QuadLaw(...
            prams1,coll_inputs,zz0,VV0,prams2)
%% need: om,d=draft,cD_dim,h,RAO_surge
%% prams1      = [om,d,RAO_surge,h]; d=draft, h=thickness
%% coll_inputs = [A,rest_coeff,cD_dim];
%% prams2      = [conc,cg,f_coll] to convert to attenuation coefficient
%% also:
%% zz0: time (\omega*t) at moment of collision
%% VV0: velocity at moment of collision
%%
%% drag law: tau=-rhow*cD*|v'|v'
%%  > cD non-dimensional
%%  > v' = v_disc-v_wave,
%%          where v_wave is the velocity the disc would have been
%%          moving at if it hadn't had a collision

do_test  = 0;
rhow     = 1025;
g        = 9.81;

%%outputs (1 for each floe);
T_relax  = [0;0];%% time till drag stress is 1st zero
WW       = [0;0];%% work done (per unit area) by drag stress
tau_init = [0;0];%% initial (max) stress 

if nargin==0
   'using test inputs'
   do_test     = 1;
   %cD          = 1e-3;%% - nondimensional
   cD          = 1e2;%% - nondimensional
   A           = 5e-2;
   h           = 3.3e-2;
   d           = 1.8e-2;%%draft
   rest_coeff  = .9;

   %% =====================================
   %%inputs to fn_attn_coeff_collisions
   coll_inputs0 = [A,rest_coeff];
   %%
   if 1
      period   = 1.1;
      RAOsurge = 0.568508129170096;
   else
      period   = 1.5;
      RAOsurge = 0.851899193743405;
   end
   om = 2*pi/period;
   k  = om^2/g;
   cg = om/2/k;
   wave_pram0 = [om,k,cg]

   sep   = 1;
   D  = .99;
   rhoi  = d/h*rhow;
   ice_pram0 = [D,h,rhoi];

   conc     = .77;
   om       = 2*pi/period;
   f_coll   = om/pi;%%2*f

   out = fn_attn_coeff_collisions(...
      coll_inputs0,wave_pram0,ice_pram0,...
      RAOsurge,conc,sep,f_coll)

   omt1     = out.omt1;
   V1       = out.V1;
   omt2     = out.omt2;
   V2       = out.V2;
   clear wave_pram0,ice_pram0,coll_inputs0;
   %% =====================================

   VV0   = [V1;V2];
   zz0   = [omt1;omt2];
   clear omt1 omt2 V1 V2;
   
   prams1      = [om,d,RAOsurge,h];
   coll_inputs = [A,rest_coeff,cD];
   %%
   prams2   = [f_coll,conc,cg];
   clear om d RAOsurge h
   clear A rest_coeff cD
   clear f_coll conc cg
end


%% prams1 = [om,d,RAO_surge,h]; d=draft
om       = prams1(1);
d        = prams1(2);
RAOsurge = prams1(3);
h        = prams1(4);

%% coll_inputs = [A,rest_coeff,cD_dim];%%NB currently scalars
A          = coll_inputs(1);
rest_coeff = coll_inputs(2);
cD         = coll_inputs(3);

%% prams2   = [f_coll,conc,cg];
f_coll = prams2(1);
conc   = prams2(2);
cg     = prams2(3);

%% nondim:
%% speed: vfac=X*omega
%% time:  z=t*omega
%%
%% drag eqn:
%% rhoi*h*dU/dt=rhoi*h*vfac*omega*dU'/dz=tau_fac*tau
%% where tau_fac=rhoi*h*X*omega^2
%% also tau=rhow*cD*vfac^2*|V'-U'|(V'-U')
%% so tau/tau_fac = (rhow*vfac^2/tau_fac)*cD*|V'-U'|(V'-U')
%%                = (X/d)*cD*|V'-U'|(V'-U')

X        = RAOsurge*A;
sd       = d*om^2/g;
cD       = cD*X/d;    %% rescaled drag coeff
Wfac     = sd*(rhow*g*X^2);   %% convert work from nondimensional to dimensional
tau_fac  = sd*(rhow*g*X);     %% convert stress from nondimensional to dimensional

for j=1:2

   V0    = VV0(j)/(X*om);%%initial speed [nondimensional]
   v0    = -rest_coeff*V0;
   z0    = zz0(j);%%initial value of z=\omega*t
   %[V0,-sin(z0)]

   %zvec  = z0+(0:N)'*dz;
   %sgn   = sign(V0);

   if 1%%use matlab solver
      zspan    = z0+[0,pi];%interval of (t\omega)
      %solver   = @ode113;col = 'k';
      %solver   = @ode23;col = 'm';
      solver   = @ode45;col   = 'b';

      %% solve rhoi*dv/dt=tau, v(0)=v0
      %% nondim u=v/(X*omega), z=t*omega
      %% du/dz=tau/tau_fac
      [zout,v1_out] = feval(solver,@(z,v1) stress_fun_quad(z,v1,cD),zspan,v0);
      V1_out        = -sin(zout);
      tau_out       = tau_fac*stress_fun_quad(zout,v1_out,cD);%%[Pa]
      %%
      tau0 = tau_out(1);
      %[tau0,tau_fac*cD*abs(V0-v0)*(V0-v0)]
      %[v0,v1_out(1)]
      %error('hey')
      tau_init(j) = tau0;

      %%find 1st time stress is 0 - defines a 'relaxation' time
      jz = find(tau0*tau_out<0,1,'first');
      if isempty(jz)
         T_relax(j) = NaN;
      else
         T_relax(j) = zout(jz)/om;
      end

      %%determine when to calculate the work done
      if 0;
         %% up to T_relax;
         jR = jz;
      else
         %% up to next collision;
         zcoll = om/f_coll;
         jR    = find((zout-z0)>zcoll,1,'first');
         if isempty(jR)
            jR = length(zout);
         end
      end
      %%
      W_out = cumtrapz(zout,work_fun_quad(zout,v1_out,cD));%%nondimensional work
      W_out = Wfac*W_out;                             %%[J/m^2]
      WW(j) = W_out(jR);
      %%
      if do_test
         figure(200+j);
         subplot(3,1,1);
         plot(zout/pi,V1_out,'r');
         hold on;
         plot(zout/pi,v1_out,col);
         yl = get(gca,'ylim');
         plot(om*T_relax(j)/pi+0*yl,yl,'g');
         hold off;
         GEN_proc_fig('\omega{t}/\pi','v_1');
         xlim(zspan/pi);
         legend('V_{wave}','V_{floe}','T_{relax}')
         %%
         subplot(3,1,2);
         tau_out  = stress_fun_quad(zout,v1_out,cD);
         plot(zout/pi,tau_out);
         xl = get(gca,'xlim');
         x  = xl(1)+.2*(xl(2)-xl(1));
         yl = get(gca,'ylim');
         y  = yl(1)+.4*(yl(2)-yl(1));
         str   = ['T_{relax}=',num2str(T_relax(j)),'s; ',...
                  '\tau_{init}=',num2str(tau_init(j)),'Pa; ',...
                  'W=',num2str(WW(j)),' Jm^{-2}'];
         txt   = text(x,y,str);
         GEN_font(txt);
         hold off;
         GEN_proc_fig('\omega{t}/\pi','\tau, Pa');
         xlim(zspan/pi);
         %%
         subplot(3,1,3);
         plot(zout/pi,W_out);
         GEN_proc_fig('\omega{t}/\pi','W, Jm^{-2}');
         xlim(zspan/pi);
         pause;
         %return;
      end
   end

   %% %% initial condition for u:
   %% %% u'(0)=-1/cD*(v0-V0)=(1-e)/cD*V0*u(0)=ic_coeff*u(0)
   %% %% u'=w => w(0)=ic_coeff*u(0)
   %% %% w'=-cos(z)/(cD*sgn)*u=f*u
   %% ic_coeff = (1-rest_coeff)/(cD*sgn)*V0;
   %% uw       = [1;ic_coeff];%%u(0) is arbitrary
   %% %%
   %% uvec     = 0*zvec;
   %% wvec     = 0*zvec;
   %% fvec     = 0*zvec;
   %% %%
   %% n           = 0;
   %% V1_vec      = -sin(zvec);
   %% v1_vec      = 0*zvec;
   %% tau_vec     = 0*zvec;
   %% tau_vec_    = 0*zvec;
   %% Wvec_       = 0*zvec;
   %% v1_vec(n+1) = v0;

   %% %% weights for trapezoidal integration (linear interpolation);
   %% wt0   = dz/2;%%weight for LH point
   %% wt1   = dz/2;%%weight for RH point

   %% W              = 0;%%work done [J/m^2] (integral of stress*v1*dt
   %% y0             = v1_vec(n+1)-V1_vec(n+1);
   %% tau            = cD*sgn*y0^2;
   %% tau_vec(n+1)   = tau;
   %% ig0            = tau*v1_vec(n+1);
   %%    %%integrand for work integral at LH of interval;

   %% z           = zvec(n+1);
   %% f           = -cos(z)/(cD*sgn);
   %% fvec(n+1)   = f;

   %% %%simple solution;
   %% U_vec          = 0*zvec;
   %% U              = v1_vec(n+1);
   %% Uvec(n+1)      = U;
   %% tau_           = tau;
   %% tau_vec_(n+1)  = tau_;
   %% ig0_           = ig0;
   %% W_             = 0;

   %% for n=1:N
   %%    M  = [0,1;f,0];
   %%    V1 = V1_vec(n+1);

   %%    %% step forward in time
   %%    uw          = uw+(M*uw)*dz;
   %%    uvec(n+1)   = uw(1);
   %%    wvec(n+1)   = uw(2);

   %%    %% alternative solution method
   %%    U              = U+tau_*dz;%%v1
   %%    % {V1-U,tau_,U,V1}
   %%    Uvec(n+1)      = U;
   %%    tau_vec_(n+1)  = tau_;
   %%    y_             = V1_vec(n+1)-U;
   %%    sgn_           = sign(y_);
   %%    tau_           = cD*sgn_*y_^2;
   %%    ig1_           = tau_*U;
   %%    W_             = Wfac*(wt0*ig0_+wt1*ig1_);
   %%    Wvec_(n+1)     = W_;

   %%    %% get v1 and check if restoring force changes direction
   %%    y1             = -cD*sgn*uw(2)/uw(1);%%v1-V1
   %%    sgn            = sign(y1);
   %%    v1             = V1+y1;
   %%    v1_vec(n+1)    = v1;
   %%    tau            = cD*sgn*y1^2;
   %%    tau_vec(n+1)   = tau;


   %%    %% Calculate work done
   %%    ig1         = tau*v1;%%integrand for work integral;
   %%    W           = W+Wfac*(wt0*ig0+wt1*ig1);
   %%    Wvec(n+1)   = W;

   %%    %%reset these for next time step:
   %%    z           = zvec(n+1);
   %%    f           = -sgn/cD*cos(z);
   %%    fvec(n+1)   = f;
   %%    ig0         = ig1;

   %%    if 1
   %%       %%plot v1 as we go:
   %%       subplot(4,1,1);
   %%       plot(zvec/pi,V1_vec,'r');
   %%       hold on;
   %%       %%
   %%       nn = 1+(0:n)';
   %%       zp = zvec(nn);
   %%       plot(zp/pi,v1_vec(nn));
   %%       plot(zp/pi,Uvec(nn),'m');
   %%       GEN_proc_fig('\omega{t}/\pi','v_1');
   %%       xlim(z0/pi+[0,2]);
   %%       hold off;

   %%       %%plot tau as we go:
   %%       subplot(4,1,2);
   %%       plot(zp/pi,tau_vec(nn));
   %%       hold on;
   %%       plot(zp/pi,tau_vec_(nn),'m');
   %%       GEN_proc_fig('\omega{t}\pi','tau');
   %%       xlim(z0/pi+[0,2]);

   %%       %%plot W as we go:
   %%       subplot(4,1,3);
   %%       plot(zp/pi,Wvec(nn));
   %%       hold on;
   %%       plot(zp/pi,Wvec_(nn),'m');
   %%       GEN_proc_fig('\omega{t}\pi','W, Jm^{-2}');
   %%       xlim(z0/pi+[0,2]);



   %%       %%plot u as we go:
   %%       subplot(4,1,4);
   %%       plot(zp/pi,uvec(nn));
   %%       GEN_proc_fig('\omega{t}\pi','u');
   %%       xlim(z0/pi+[0,2]);
   %%       %pause(.1);
   %%    end
   %% end
end
alp = conc*f_coll*sum(WW)/cg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tau   = stress_fun_quad(t,v1,cD)
%% stress from drag law
%% - non-dimensional

V1  = -sin(t);
y   = v1-V1;
sgn = -sign(y);
tau = cD*sgn.*y.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W   = work_fun_quad(t,v1,cD)
%% integrand stress*[dx/dt]
%% since \int stress*dx=\int stress*[dx/dt]*dt
%% dx/dt=relative velocity, y=v1-V1

V1    = -sin(t);
y     = v1-V1;
tau   = stress_fun_quad(t,v1,cD);
%W     = abs(tau.*v1);
W     = abs(tau.*y);%%relative velocity
