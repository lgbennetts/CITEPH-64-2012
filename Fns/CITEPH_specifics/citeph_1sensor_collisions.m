%function [an_mat,Sn_mat,t_int] = citeph_1sensor_collisions(time,displ,T_target,H_target,outloc,inloc)
function citeph_1sensor_collisions(time,displ,T_target,H_target,outloc,inloc)
%% citeph_1sensor_movingFFT.m
%% Author: Timothy Williams
%% Date:   20130722, 09:20:07 CEST

DO_PLOT     = 1;
DO_SAVE     = 0;
data_type   = 1;%%default (1) is wave elevation, 2->acceleration
%%
expt_name      = '';
sensor_name    = '';


if nargin==0
   basedir  = citeph_user_specifics;
   fdir     = [basedir '/results_preliminary/conc_79/regular/'];
   fname    = [fdir,'24071018.a13/accelerometers/houle_reg_060.dat']%A4x
   ylab     = 'a_x, ms^{-2}';
   A        = load(fname);{A}
   %%
   scale_fac   = 100;
   time        = A(:,1)/sqrt(scale_fac);
   displ       = A(:,2)/scale_fac;
   %%
   T_target = .65;
   DO_PLOT  = 1;
   %%
   if 0
      outloc   = {[basedir '/results_preliminary/conc_79/regular/23071525.a13_processed/'],...
                    '23071525.a13',...
                    '/Az/',...
                    'A1z'};
   end
end

if exist('outloc')
   DO_SAVE  = 1;
   %%
   dir1  = outloc{1}
   if ~exist(dir1)
      eval(['!mkdir ' dir1]);
   end
   %%
   expt_name   = outloc{2};%%date,time,year;
   %outdir0     = dir1;
   outdir0     = [dir1 '/' expt_name]
   if ~exist(outdir0)
      eval(['!mkdir ' outdir0]);
   end
   %%
   dir3     = outloc{3};%% eg Az
   outdir   = [outdir0 '/' dir3]
   if ~exist(outdir)
      eval(['!mkdir ' outdir]);
   end
   %%
   sensor_name = outloc{4};%%eg S1,...,S20, A1x, A1y, A1z,...,A6x, A6y, A6z
   outdir2  = [outdir,'/mat_files/'];
   if ~exist(outdir2)
      eval(['!mkdir ' outdir2]);
   end
   outfile     = [outdir2,'/collisions_',sensor_name,'.mat'];
   %%
   outdir2  = [outdir,'/eps_files/'];
   if ~exist(outdir2)
      eval(['!mkdir ' outdir2]);
   end
   figfile1    = [outdir2,'/collisions_',sensor_name,'.eps'];
   figfile2    = [outdir2,'/collisions_hist_',sensor_name,'.eps'];
   figfile3    = [outdir2,'/collisions_spec_',sensor_name,'.eps'];
   %%
   if strcmp(sensor_name(3),'x')
      ylab  = 'a_x, ms^{-2}';
   elseif strcmp(sensor_name(3),'y')
      ylab  = 'a_y, ms^{-2}';
   elseif strcmp(sensor_name(3),'z')
      ylab  = 'a_z, ms^{-2}';
   end
   %figfile2    = [outdir,'/collisions_',sensor_name,'.mat'];

   if strcmp(sensor_name(1),'A')
      data_type   = 2;
   end
end


N     = length(time);
%%
dt    = time(2)-time(1);
T     = time(end)-time(1);
fs    = 1/dt;
fmax  = fs/2;
%figure(1),plot(time,disp);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart      = 50;%%wave maker starts [s]
dist        = 13.244;%%dist to wave maker [m]
om          = 2*pi/T_target;
waveno      = om^2/9.81;
groupvel    = om/waveno/2;%%group velocity [m/s] is half the phase velocity at infinite depth
trav_time   = dist/groupvel;
ts          = tstart+trav_time
%%
tstop = T_target*120;%%approx time wave maker stops
tf    = tstop+trav_time;

j_rel = find(time>=ts & time<=tf);%%when wave-maker starts
if mod(length(j_rel),2)==1
   j_rel(1) = [];
end
t_rel    = time(j_rel);
a_rel0   = displ(j_rel);
a_rel    = a_rel0-mean(a_rel0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1%%plot histogram
   figure(101);
   [N,X] = hist(a_rel,100);
   xlab2 = 'a, ms^{-2}';
   dx    = X(2)-X(1);
   A     = dx*sum(N);
   Y     = N/A;%%normalise so area=1;
   %%
   [N2,X2]  = hist((a_rel).^2,100);
   xlab2b   = 'a^2, m^2s^{-4}';
   dx2      = X(2)-X(1);
   A2       = dx*sum(N2);
   Y2       = N2/A2;%%normalise so area=1;
   %%
   nr = 2;
   subplot(nr,2,1);
   plot(X,Y);
   GEN_proc_fig(xlab2,'Probability density');
   ttl   = title(['Acceleration histograms (',sensor_name,')']);
   GEN_font(ttl);

   subplot(nr,2,2);
   yy = dx*fliplr(cumsum(fliplr(Y)));
   plot(X,yy);
   set(gca,'yscale','log','xscale','log');
   GEN_proc_fig(xlab2,'Cumulative probability');

   if nr==2
      subplot(nr,2,3);
      plot(X2,Y2);
      GEN_proc_fig(xlab2b,'Probability density');

      subplot(nr,2,4);
      yy = dx2*fliplr(cumsum(fliplr(Y2)));
      plot(X2,yy);
      set(gca,'yscale','log','xscale','log');
      GEN_proc_fig(xlab2b,'Cumulative probability');
   end

   saveas(gcf,figfile2,'epsc');
end


T_filt   = T_target/4%%s
f_filt   = 1/T_filt;%%Hz
[Bh,Ah]  = butter(3,f_filt/fmax,'high');
[Bl,Al]  = butter(3,f_filt/fmax,'low');
%%
ah = filter(Bh,Ah,a_rel);
al = filter(Bl,Al,a_rel);
%%
[Sn_a ,ff,an_cx,t_mean,xtra]  = citeph_get_ft(t_rel,a_rel);T_target,xtra
asig1 = xtra{1}/2;
Tp1   = xtra{2}
[Sn_ah,ff,an_cx,t_mean,xtra]  = citeph_get_ft(t_rel,ah);
asig2    = xtra{1}/2
Tp2   = xtra{2}
[Sn_al,ff,an_cx,t_mean,xtra]  = citeph_get_ft(t_rel,al);
asig3 = xtra{1}/2;
Tp3   = xtra{2}

a_thresh = 2*asig2;%% 4 std deviations in amp of acc
jcol     = find(abs(ah)>a_thresh);
tcol     = t_rel(jcol);
acol     = ah(jcol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1%%if time between big accelerations is too small, they are the same collision
   Tc    = 10*dt;%%max collision time
   %%
   tcol2a         = [];
   tcol2b         = [];
   tcol_current   = [];
   acol2          = [];
   for j=1:length(tcol)
      t0 = tcol(j);
      t1 = t0+Tc;
      jj = find(t0<tcol & tcol<=t1);

      tcol_current   = [tcol_current;t0];
      if isempty(jj)%%collision finished
         t_a            = tcol_current(1);
         t_b            = t0;
         tcol2a         = [tcol2a;t_a];
         tcol2b         = [tcol2b;t_b];
         %%
         jc    = find(t_a<=tcol & tcol<=t_b);
         acol2 = [acol2;sqrt( mean(acol(jc).^2) )];
         %%
         tcol_current   = [];
      end
      %tcol_current,pause
   end

   tcol2    = (tcol2a+tcol2b)/2;
   col_dt   = tcol2a-tcol2b;
else
   tcol2    = tcol;
   col_dt   = dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_SAVE==1
   save(outfile,'tcol2','acol2','col_dt','asig1','asig2','asig3','T_target','H_target');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103);
subplot(3,1,1), plot(t_rel,a_rel);
j0       = strfind(inloc,'conc');
inloc0   = inloc(j0:end);
inloc0   = strrep(inloc0,'_','\_');

ttl   = title({['Raw time series (',sensor_name, '): a_s = ',num2str(asig1,'%0.2e'),' ms^{-2}; ',...
               'T_p = ',num2str(T_target), ' s'],inloc0});
GEN_font(ttl);
if 0
   yl = .3*[-1 1];
   ylim(yl);
else
   yl = get(gca,'ylim');
end
GEN_proc_fig('Time, s',ylab)
%%
subplot(3,1,3), plot(t_rel,al);
ttl   = title(['Low-pass-filtered time series: ',num2str(100*(asig3/asig1)^2,'%0.1f'),' % of variance']);
GEN_font(ttl);
GEN_proc_fig('Time, s',ylab)
%%
subplot(3,1,2);
plot(t_rel,ah);%%plot once 1st to get y-range
if 0
   ylim(yl);
else
   yl = get(gca,'ylim');
end
xl = get(gca,'xlim');
%%
for n=1:length(tcol)
   plot(tcol(n)+0*yl,yl,'g');
   if n==1
      hold on;
   end
end
hold on;
plot(xl,0*xl+a_thresh,'m');
plot(xl,0*xl-a_thresh,'m');
plot(t_rel,ah);%%plot again so it's on top of other lines
hold off;
ttl   = title(['High-pass-filtered time series: ',num2str(100*(asig2/asig1)^2,'%0.1f'),' % of variance']);
GEN_font(ttl);
GEN_proc_fig('Time, s',ylab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
   saveas(gcf,figfile1,'epsc');
end

if 1%%plot spectra
   figure(104);
   subplot(3,1,1), plot(1./ff,Sn_a);
   set(gca,'yscale','log');
   hold on;
   yl = get(gca,'ylim');
   plot(T_target+0*yl,yl,'r');
   plot(Tp1+0*yl,yl,'--k');
   plot(T_filt+0*yl,yl,'g');
   xlim([0 15]);
   GEN_proc_fig('Period, s','SD, m^2s^{-3}');
   ttl   = title(['Spectral density (',sensor_name,')']);
   GEN_font(ttl);
   hold off;
   %%
   subplot(3,1,2), plot(1./ff,Sn_ah);
   [Smax,jmax]  = max(Sn_ah)
   set(gca,'yscale','log');
   hold on;
   yl = get(gca,'ylim')
   plot(T_target+0*yl,yl,'r');
   plot(Tp2+0*yl,yl,'--k');
   plot(T_filt+0*yl,yl,'g');
   xlim([0 2*T_target]);
   xl = get(gca,'xlim')
   plot(xl,0*xl+Smax,'k')
   hold off;
   GEN_proc_fig('Period, s','SD (HPF), m^2s^{-3}');
   %%
   subplot(3,1,3), plot(1./ff,Sn_al);
   set(gca,'yscale','log');
   hold on;
   yl = get(gca,'ylim');
   plot(T_target+0*yl,yl,'r');
   plot(Tp3+0*yl,yl,'--k');
   plot(T_filt+0*yl,yl,'g');
   xlim([0 15]);
   hold off;
   GEN_proc_fig('Period, s','SD (LPF), m^2s^{-3}');

   if 1
      saveas(gcf,figfile3,'epsc');
   end
end
