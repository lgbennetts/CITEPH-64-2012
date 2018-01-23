function [alp,WW,tau_init,T_relax] = fn_ode_drag_LinearLaw(...
            prams1,coll_inputs,zz0,VV0,prams2)
%% need: om,d=draft,cD_dim,h,RAO_surge
%% prams1      = [om,d,RAO_surge,h]; d=draft, h=thickness
%% coll_inputs = [A,rest_coeff,cD];
%% prams2      = [conc,cg,f_coll] to convert to attenuation coefficient
%% also:
%% zz0: time (\omega*t) at moment of collision
%% VV0: velocity at moment of collision
%%
%% drag law: tau=-(rhow*d*omega)*kD*v'
%%  > kD non-dimensional
%%  > v' = v_disc-v_wave,
%%          where v_wave is the velocity the disc would have been
%%          moving at if it hadn't had a collision
%%  > calibrate kD so initial stress is the 
%%       same as initial stress from quadratic drag law:
%%       tau=-rhow*cD*|v'|v'
%%  > ie kD=cD*|v'|/(d*omega)
%% pass in cD

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
A           = coll_inputs(1);
rest_coeff  = coll_inputs(2);%%non-dimensional
cD          = coll_inputs(3);%%non-dimensional drag coeff - scaled by \rho_water*d*\omega
%coll_inputs,pause
%% prams2   = [f_coll,conc,cg];

f_coll   = prams2(1);
conc     = prams2(2);
cg       = prams2(3);

X        = RAOsurge*A;
sd       = om^2/g*d;
cD       = cD*X/d %%rescaled drag coeff: du/dz=cD*|V-v|*(V-v)
Wfac     = sd*(rhow*g*X^2);  %%convert work from nondimensional to dimensional
E        = rhow*g*A^2/2;
Wfac0    = 2*RAOsurge^2*sd;%[Wfac0*E,Wfac],pause 
tau_fac  = sd*(rhow*g*X);    %%convert stress from nondimensional to dimensional

for j=1:2

   V0    = VV0(j)/(X*om);%%initial speed [nondimensional]
   v0    = -rest_coeff*V0;
   z0    = zz0(j);%%initial value of z=\omega*t
   zspan = z0+[0,om/f_coll];

   %% calibrate linear drag law
   %% using cD and initial speed
   %% (this gives the same initial stress
   %%  as from a quadratic drag law)
   kD = cD*abs(V0-v0);
   %%
   eps         = 5e-4;
   z_relax     = -log(eps)/kD;
   T_relax(j)  = z_relax/om;
   %%
   N     = 2e2;

   zvec  = z0+linspace(0,z_relax,N)';%%more points at start
   zzz   = linspace(zvec(end),zspan(2),N)';
   zvec  = [zvec;zzz(2:end)];
   clear zzz;

   f1    = particular_soln(zvec,kD);
   Cj    = v0-f1(1);
   u1    = Cj*exp(-kD*(zvec-z0))+f1;
   tau1  = tau_fac*stress_fun_lin(zvec,u1,kD);%%dimensional stress
   w1    = cumtrapz(zvec,work_fun_lin(zvec,u1,kD));%%nondimensional work
   W1    = Wfac*w1;%%dimensional work
   WW(j) = Wfac0*trapz(zvec,work_fun_lin(zvec,u1,kD));%%work/E->attenuation coefficient

   tau0 = tau1(1)
   [v0,u1(1)]
   tau_init(j) = tau0;

   if do_test%%plot and compare with ode45
      tst_work = [WW(j)*E,W1(end)]

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      figure(300+j);
      subplot(3,1,1);
      plot(zvec/pi,u1,'k');%%analytic soln
      hold on;
      %plot(zvec/pi,f1,'r--');
      xlim(zspan/pi);
      subplot(3,1,2);
      plot(zvec/pi,tau1,'k');
      hold on;
      xlim(zspan/pi);
      subplot(3,1,3);
      plot(zvec/pi,W1,'k');
      hold on;
      xlim(zspan/pi);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %solver   = @ode113;col = 'k';
      %solver   = @ode23;col = 'm';
      solver   = @ode45,col   = '--r';
      %%
      [zout,v1_out]  = feval(solver,@(z,v1) stress_fun_lin(z,v1,kD),zspan,v0);
      V1_out         = -sin(zout);
      tau_out        = tau_fac*stress_fun_lin(zout,v1_out,kD);%%[Pa]
      if 0
         figure(4)
         plot(zout,tau_fac*stress_fun_lin(zout,v1_out,kD),'b',...
              zvec,tau_fac*stress_fun_lin(zvec,u1,kD),'k');
         {zvec(1),zout(1)}
         pause
         close
      end
      %%

      % %%find 1st time stress is 0 - defines a 'relaxation' time
      % jz          = find(tau0*tau_out<0,1,'first');
      % T_relax(j)  = zout(jz)/om;

      % %%determine when to calculate the work done
      % if 0;
      %    %% up to T_relax;
      %    jR = jz;
      % else
      %    %% up to next collision;
      %    zcoll = om/f_coll;
      %    jR    = find((zout-z0)>zcoll,1,'first');
      %    if isempty(jR)
      %       jR = length(zout);
      %    end
      % end
      % %%
      W_out = cumtrapz(zout,work_fun_lin(zout,v1_out,kD));%%nondimensional work
      W_out = Wfac*W_out;                             %%[J/m^2]
      % WW(j) = W_out(jR);
      if 0
         figure(4)
         ig1   = Wfac*work_fun_lin(zvec,u1,kD);
         ig2   = Wfac*work_fun_lin(zout,v1_out,kD);{zvec,ig1;zout,ig2}
         plot(zvec,ig1,'b',...
              zout,ig2,'k');
         {zvec(1),zout(1);
         trapz(zvec,ig1),trapz(zout,ig2);
         trapz(zvec,W1),trapz(zout,W_out)}
         pause
         close
      end
      %%
      figure(300+j);
      subplot(3,1,1);
      plot(zout/pi,V1_out,'b');
      hold on;
      plot(zout/pi,v1_out,col);%%ode45 soln
      yl = get(gca,'ylim');
      plot((z0+z_relax)/pi+0*yl,yl,'g');
      hold off;
      GEN_proc_fig('\omega{t}/\pi','v_1');
      xlim(zspan/pi);
      legend('V_{floe} (AS)','V_{wave}','V_{floe} (NS)','T_{relax}')
      %%
      subplot(3,1,2);
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
   end%%end do_test
end%%end j loop

alp   = conc*f_coll*sum(WW)/cg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=particular_soln(z,kD)

f  = -kD^2/(kD^2+1)*sin(z) + kD/(kD^2+1)*cos(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tau   = stress_fun_lin(t,u1,kD)
%%stress from drag law

v1    = -sin(t);
y     = u1-v1;
tau   = -kD*y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W   = work_fun_lin(t,u1,kD)
%% integrand stress*[dx/dt]
%% since \int stress*dx=\int stress*[dx/dt]*dt
%% dx/dt=relative velocity, y=v1-V1

v1    = -sin(t);
y     = u1-v1;
tau   = stress_fun_lin(t,u1,kD);
W     = abs(tau.*y);%%relative velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tau   = stress_fun_quad(t,v1,cD)
%% stress from quadratic drag law
%% - use this to calibrate linear drag coefficient

V1    = -sin(t);
y     = v1-V1;
sgn   = -sign(y);
tau   = cD*sgn.*y.^2;
