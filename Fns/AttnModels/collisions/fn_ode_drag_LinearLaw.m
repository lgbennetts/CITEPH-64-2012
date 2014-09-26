function [alp,WW,tau_init,T_relax] = fn_ode_drag_LinearLaw(...
            prams1,coll_inputs,zz0,VV0,prams2)
%% need: om,d=draft,cD_dim,h,RAO_surge
%% prams1      = [om,d,RAO_surge,h]; d=draft
%% coll_inputs = [A,rest_coeff,cD_dim];
%% prams2      = [conc,cg,f_coll] to convert to attenuation coefficient
%% also:
%% zz0: time (\omega*t) at moment of collision
%% VV0: velocity at moment of collision

do_test  = 0;
rhow     = 1025;
g        = 9.81;

%%outputs (1 for each floe);
T_relax  = [0;0];%% time till drag stress is 1st zero
WW       = [0;0];%% work done (per unit area) by drag stress
tau_init = [0;0];%% initial (max) stress 

if nargin==0
   'test inputs'
   do_test     = 1;
   cD          = 1e1;%%kg/m^3
   A           = 2e-2;
   h           = 3.3e-2;
   d           = 1.8e-2;%%draft
   rest_coeff  = .9;

   if 1
      om       = 2*pi/1.1;
      omt1     = 1.205531059037677;
      V1       = -0.060661652608385;
      omt2     = -1.021253028778352;
      V2       = 0.055383740125916;
      RAOsurge = 0.568508129170096;
   else
      om       = 2*pi/1.5;
      omt1     = 0.508439575729424;
      V1       = -0.034743261608009;
      omt2     = -1.861452883782580;
      V2       = 0.068375052575801;
      RAOsurge = 0.851899193743405;
   end

   VV0   = [V1;V2];
   zz0   = [omt1;omt2];

   tau_init       = 2e4;%%Pa - enable for comparison with quadratic drag law
   y_init         = V1+rest_coeff*V1;
   cD             = abs( tau_init/(om*d*rhow*y_init) )
   clear tau_init y_init;
   
   prams1      = [om,d,RAOsurge,h];
   coll_inputs = [A,rest_coeff,cD];
   %%
   f_coll   = om/pi;%%2*f
   conc     = .77;
   cg       = 1;
   prams2   = [f_coll,conc,cg];
   clear om d RAOsurge h
   clear A rest_coeff cD
   clear f_coll conc cg
   clear omt1 omt2 V1 V2;
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
Wfac     = sd*(rhow*g*X^2);  %%convert work from nondimensional to dimensional
E        = rhow*g*A^2/2;
Wfac0    = 2*RAOsurge^2*sd;%[Wfac0*E,Wfac],pause 
tau_fac  = sd*(rhow*g*X);    %%convert stress from nondimensional to dimensional

for j=1:2

   eps         = 5e-4;
   z_relax     = -log(eps)/cD;
   T_relax(j)  = z_relax/om;
   %%
   V0    = VV0(j)/(X*om);%%initial speed [nondimensional]
   v0    = -rest_coeff*V0;
   z0    = zz0(j);%%initial value of z=\omega*t
   zspan = z0+[0,om/f_coll];
   %%
   N     = 2e2;

   zvec  = z0+linspace(0,z_relax,N)';%%more points at start
   zzz   = linspace(zvec(end),zspan(2),N)';
   zvec  = [zvec;zzz(2:end)];
   clear zzz;

   f1    = particular_soln(zvec,cD);
   Cj    = v0-f1(1);
   u1    = Cj*exp(-cD*(zvec-z0))+f1;
   tau1  = tau_fac*stress_fun(zvec,u1,cD);%%dimensional stress
   w1    = cumtrapz(zvec,work_fun(zvec,u1,cD));%%nondimensional work
   W1    = Wfac*w1;%%dimensional work
   WW(j) = Wfac0*trapz(zvec,work_fun(zvec,u1,cD));%%work/E->attenuation coefficient

   if do_test%%plot and compare with ode45
      tst_work = [WW(j)*E,W1(end)]

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      figure(2);
      subplot(3,1,1);
      plot(zvec/pi,u1,'k',zvec/pi,-sin(zvec),'r');
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
      pause
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %solver   = @ode113;col = 'k';
      %solver   = @ode23;col = 'm';
      solver   = @ode45,col   = '--b';
      %%
      [zout,v1_out]  = feval(solver,@(z,v1) stress_fun(z,v1,cD),zspan,v0);
      V1_out         = -sin(zout);
      tau_out        = tau_fac*stress_fun(zout,v1_out,cD);%%[Pa]
      if 0
         figure(4)
         plot(zout,tau_fac*stress_fun(zout,v1_out,cD),'b',...
              zvec,tau_fac*stress_fun(zvec,u1,cD),'k');
         {zvec(1),zout(1)}
         pause
         close
      end
      %%
      tau0        = tau_out(1);
      tau_init(j) = tau0;

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
      W_out = cumtrapz(zout,work_fun(zout,v1_out,cD));%%nondimensional work
      W_out = Wfac*W_out;                             %%[J/m^2]
      % WW(j) = W_out(jR);
      if 0
         figure(4)
         ig1   = Wfac*work_fun(zvec,u1,cD);
         ig2   = Wfac*work_fun(zout,v1_out,cD);{zvec,ig1;zout,ig2}
         plot(zvec,ig1,'b',...
              zout,ig2,'k');
         {zvec(1),zout(1);
         trapz(zvec,ig1),trapz(zout,ig2);
         trapz(zvec,W1),trapz(zout,W_out)}
         pause
         close
      end
      %%
      figure(2);
      subplot(3,1,1);
      plot(zout/pi,V1_out,'r');
      hold on;
      plot(zout/pi,v1_out,col);
      yl = get(gca,'ylim');
      plot((z0+z_relax)/pi+0*yl,yl,'g');
      hold off;
      GEN_proc_fig('\omega{t}/\pi','v_1');
      xlim(zspan/pi);
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
function f=particular_soln(z,cD)

f  = -cD^2/(cD^2+1)*sin(z) + cD/(cD^2+1)*cos(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tau   = stress_fun(t,u1,cD)
%%stress from drag law

v1    = -sin(t);
y     = u1-v1;
tau   = -cD*y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W   = work_fun(t,u1,cD)
%% integrand stress*[dx/dt]
%% since \int stress*dx=\int stress*[dx/dt]*dt
%% dx/dt=relative velocity, y=v1-V1

v1    = -sin(t);
y     = u1-v1;
tau   = stress_fun(t,u1,cD);
W     = abs(tau.*y);%%relative velocity
