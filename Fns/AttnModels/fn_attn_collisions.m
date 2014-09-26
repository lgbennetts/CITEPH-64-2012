%% attn_collisions.m
%% Author: Timothy Williams
%% Date: 20140910, 15:40:45 CEST

function out   =...
   fn_attn_collisions(coll_inputs,wave_pram0,ice_pram,tank_pram,scat_pram,conc,sep)

g        = 9.81;
rhow     = 1025;
do_test  = 0;
SAVE_FIG = 1;%%save figure (applies only if do_test==1)

if nargin==0
   do_test  = 1;
   %%
   if 0
      om       = 2*pi/.6;
      alp_scat = 1.338103972350944;
      RAOsurge = 0.100165667108768;
   elseif 1
      om       = 2*pi/.7;
      alp_scat = 0.747557856183909;
      RAOsurge = 0.052160250017141;
   elseif 1
      om       = 2*pi/1.1;
      alp_scat = 0.010949152443493;
      RAOsurge = 0.568508129170096;
   elseif 1
      om       = 2*pi/1.5;
      alp_scat = 3.766326568716634e-05;
      RAOsurge = 0.851899193743405;
   end
   k           = om^2/g;
   cp          = om/k;
   wave_pram0  = [om,k];
   clear om k
   %%
   if 0
      A0          = 50e-2;%+zeros(5,1);
      rest_coeff  = .9+0*A0;
      coll_inputs = [A0,rest_coeff];
   elseif 0
      coll_inputs =...
        [10e-2    0.9
         15e-2    0.9
         20e-2    0.9
         40e-2    0.9
         50e-2    0.9];
   else
      %%add drag
      A0          = 20e-2;%+zeros(5,1);
      rest_coeff  = .9+0*A0;
      cD          = 1e4+0*A0;%%non-dim (linear law)
      coll_inputs = [A0,rest_coeff,cD];
      clear cD;
   end
   clear A0 rest_coeff
   %%
   D           = .99;
   sep         = 1.00;
   conc        = (pi*D^2/4)/(sep^2);
   h           = 33e-3;
   rhoi        = 18/33*rhow;
   ice_pram    = [D,h,rhoi];
   clear D h rhoi
   %%
   scat_pram   = [alp_scat,RAOsurge];
   clear alp_scat RAOsurge
   %%
   MIZwidth    = 5;
   H           = 3.1;
   tank_pram   = [MIZwidth,H];
   clear MIZwidth H
end
%coll_inputs,wave_pram0,ice_pram,tank_pram,scat_pram,conc,sep

MIZwidth = tank_pram(1);
H        = tank_pram(2);
%%
A0          = coll_inputs(:,1);
rest_coeff  = coll_inputs(:,2);
Ncoll       = length(A0);
DO_DRAG     = 0;
if size(coll_inputs,2)==3
   DO_DRAG  = 1;%%include drag; 3rd column is drag coefficient
   drag_fxn = @fn_ode_drag_LinearLaw; %%linear drag law
   %drag_fxn = @fn_ode_drag_QuadLaw;  %%quadratic drag law
end
%%
om       = wave_pram0(1);
k        = wave_pram0(2);
f_coll   = om/pi;%%2*f
%%
cg          = om/2/k;%%inf depth
cg          = cg+H/(2*om*g*k)*(g^2*k^2-om^4);%%finite depth correction
wave_pram   = [wave_pram0,cg];

nx = 200;
x  = linspace(0, MIZwidth,nx+1)';
dx = x(2);

alp_scat = scat_pram(1);
RAOsurge = scat_pram(2);
[alp_coll,out2]   = fn_attn_coeff_collisions(...
   coll_inputs,wave_pram,ice_pram,RAOsurge,conc,sep,f_coll);

%A_all    = A0*eye(nx+1,1);
%ac_all   = alp_coll*eye(nx+1,1);
A_all       = zeros(Ncoll,nx+1);
A_all(:,1)  = A0;
%%
ac_all      = zeros(Ncoll,nx+1);
ac_all(:,1) = alp_coll;

%%finite difference
S0    = A0.^2/2;
S     = S0;
alp   = alp_scat+alp_coll;
COLL  = ~isempty(find(alp_coll));

if DO_DRAG & COLL
   ad_all   = zeros(Ncoll,nx+1);
   alp_drag = 0*A0;
   h        = ice_pram(2);
   d        = ice_pram(3)*h/rhow;

   for jd   = 1:Ncoll
      ci0   = coll_inputs(jd,:);
      omt0  = [out2(1).value,out2(2).value];%%collision times
      vel0  = [out2(1).value,out2(2).value];%%collision velocities
      %%
      drag_prams1 = [om,d,RAOsurge,h];
      drag_prams2 = [f_coll,conc,cg];
      %%
      [ad,WW,tau_init,T_relax] = feval(drag_fxn,...
         drag_prams1,ci0,omt0,vel0,drag_prams2);
      %%
      alp_drag(jd) = ad;
   end
   ad_all(:,1) = alp_drag;
   alp         = alp+alp_drag;
end

for n=2:nx+1
   S                 = S.*exp(-alp*dx);
   A                 = sqrt(2*S);
   coll_inputs(:,1)  = A;
   %%
   if COLL
      [alp_coll,out2]   = fn_attn_coeff_collisions(...
         coll_inputs,wave_pram,ice_pram,RAOsurge,conc,sep,f_coll);
      COLL  = ~isempty(find(alp_coll));
   else
      alp_coll = 0*A0;
   end
   %%
   alp           = alp_scat+alp_coll;
   A_all(:,n)    = A;
   ac_all(:,n)   = alp_coll;

   if DO_DRAG & COLL
      alp_drag = 0*A0;
      for jd   = 1:Ncoll
         ci1   = coll_inputs(jd,:);
         omt0  = [out2(1).value,out2(2).value];%%collision times
         vel0  = [out2(1).value,out2(2).value];%%collision velocities
         %%
         [ad,WW,tau_init,T_relax] = feval(drag_fxn,...
            drag_prams1,ci1,omt0,vel0,drag_prams2);
         %%
         alp_drag(jd)   = ad;
         % {ad,alp_coll}
      end
      ad_all(:,n) = alp_drag;
      alp         = alp+alp_drag;
   end
end

%%output transmitted energy:
E_t   = real(S./S0);%%transmitted energy
out_str  = [ ' ''transmitted energy'' '];
out_val  = [ 'E_t' ];

%%output transmitted energy with no collisions:
E_t0  = real(exp(-alp_scat*MIZwidth));%%transmitted energy without collisions
out_str  = [out_str '; ''transmitted energy ratio (no coll)'' '];
out_val  = [out_val '; E_t/E_t0' ];

%%output phase factor cos(ks)
phase_fac   = cos(k*sep);%%transmitted energy
out_str     = [out_str '; ''phase factor cos(ks)'' '];
out_val     = [out_val '; phase_fac' ];

%%output wave period
out_str  = [out_str '; ''wave period'' '];
out_val  = [out_val '; 2*pi/om' ];

cmd   = ['out=struct( ''name'', {' out_str ...
         '}, ''value'', {' out_val '});'];
eval(cmd);

if 0%do_test==0
   for lj=1:length(out)
      disp([out(lj).name,':']);
      disp(out(lj).value);
   end
end

if do_test
   for lj=1:length(out)
      disp([out(lj).name,':']);
      disp(out(lj).value);
   end
   Efac  = rhow*g;
   
   %%plot amplitude;
   dd = diag(1./A0);
   subplot(2,1,1);
   plot(x,(dd*A_all).^2);
   hold on;
   plot(x,exp(-alp_scat*x),'--');
   [Am,jm]  = max(A0);
   ttl   = title(['Period: ',num2str(2*pi/om),...
                  's; Rest. Coeff.: ', num2str(coll_inputs(1,2))]);
   GEN_font(ttl);
   hold off;
   GEN_proc_fig('x, m','E/E0');
   leg   = [];
   for loop_j=1:Ncoll
      A     = A0(loop_j)*100;%%amp in cm
      rc    = coll_inputs(loop_j,2);
      leg   = [leg, ', ''(',num2str(A),'cm,',num2str(rc),')'' '];
   end
   cmd  = ['Lg = legend(',leg(2:end),...
           ',''No Collisions'', ''Location'' ,',' ''EastOutside'' );'];
   eval(cmd);

   %%plot alp_coll;
   subplot(2,1,2);
   leg   = [];
   if DO_DRAG
      drag_fxn
      plot(x,ac_all+alp_scat+ad_all,'-.');
      hold on;
      leg   = [leg, ', ''drag'' '];
   end
   plot(x,ac_all+alp_scat);
   hold on;
   plot(x,0*x+alp_scat,'--');
   GEN_proc_fig('x, m','\alpha_{tot}, m^{-1}');
   for loop_j=1:Ncoll
      A     = A0(loop_j)*100;%%amp in cm
      rc    = coll_inputs(loop_j,2);
      leg   = [leg, ', ''(',num2str(A),'cm,',num2str(rc),')'' '];
   end
   cmd  = ['Lg = legend(',leg(2:end),...
           ',''No Collisions'',  ''Location'' ,',' ''EastOutside'' );'];
   eval(cmd);

   if SAVE_FIG%%save figure
      figname  = ['eg_coll_T',num2str(2*pi/om),...
                  's_rc', num2str(rc(1)),'.eps']
      saveas(gcf,['out/',figname],'epsc');
   end
end
