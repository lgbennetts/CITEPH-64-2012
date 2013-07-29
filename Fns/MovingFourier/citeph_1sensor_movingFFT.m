%% citeph_1sensor_movingFFT.m
%% Author: Timothy Williams
%% Date:   20130722, 09:20:07 CEST
function [an_mat,Sn_mat,t_int] = citeph_1sensor_movingFFT(time,displ,T_target,outloc,inloc)

DO_PLOT     = 0;
DO_SAVE     = 0;
data_type   = 1;%%default (1) is wave elevation, 2->acceleration


if nargin==0
   fdir  = '/work/timill/CITEPH-data/results_preliminary/conc_79/regular/';
   fname = [fdir,'23071525.a13/houle_reg_061.dat']
   A     = load(fname);{A}
   %%
   scale_fac   = 100;
   time        = A(:,1)/sqrt(scale_fac);
   displ       = A(:,2)/scale_fac;
   %%
   T_target = .65;
   DO_PLOT  = 1;
   %%
   outloc   = {'/work/timill/CITEPH-data/results_preliminary/conc_79/regular/',...
                 '23071525.a13',...
                 '/Az/',...
                 'A1z'};

elseif 0
   fdir  = '/work/timill/CITEPH-data/calibration_waves/regular/';
   fname = [fdir,'18070828.a13/calib_houle_reg_061.dat']
   A     = load(fname);
   %%
   scale_fac   = 100;
   time        = A(:,1)/sqrt(scale_fac);
   displ       = A(:,2)/scale_fac;
   %%
   T_target = 2;
   DO_PLOT  = 1;
end

if exist('outloc')
   DO_SAVE     = 1;
   %%
   dir1        = outloc{1};
   if ~exist(dir1)
      eval(['!mkdir ' dir1]);
   end
   %%
   expt_name   = outloc{2};%%date,time,year;
   outdir0     = [dir1 '/' expt_name];
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
   outfile     = [outdir,'/movingFFT_',sensor_name,'.mat'];

   if strcmp(sensor_name(1),'A')
      data_type   = 2;
      %data_type   = 1;%%was getting strange results with 2
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

T_window = 10*T_target
df       = 1/T_window;
Nw       = 2^ceil(log2(T_window/dt))
T_window = dt*Nw;
t_window = (0:dt:(Nw-1))'*dt;
%%
Nf_S  = fmax*T_window
%%
t_req    = 0.5;%s
n_shift  = floor(t_req/dt);
t_req    = n_shift*dt;
%%
Nwindows = floor((T-T_window)/t_req);
Sn_mat   = zeros(Nf_S,Nwindows);
an_mat   = zeros(Nw,Nwindows);

if 0
   %m0 = floor(100/T_window)
   m0 = 10
else
   m0 = 1;
end
t_int    = zeros(Nwindows,2);
tm_vec   = zeros(Nwindows);
Tp_vec   = zeros(Nwindows);
Hs_vec   = zeros(Nwindows);

for m=m0:Nwindows
   jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
   time0       = time(jj);
   disp0       = displ(jj);
   t_start     = time0(1);
   t_int(m,:)  = time0([1 end]);
   %{time0,disp0}
   [S0,ff,an0,tm_vec(m),xtra] = get_ft(time0,disp0,data_type);
   Hs_vec(m)   = xtra{1};
   Tp_vec(m)   = xtra{2};
   Tm02_vec(m) = xtra{3};
   %{an_mat,an0}
   %{Sn_mat,S0}
   an_mat(:,m) = an0;
   Sn_mat(:,m) = S0;
   %%
   %[m,Nwindows]
   %pause;
end

if DO_SAVE
   df = ff(2)-ff(1);
   dt = time(2)-time(1);
   %disp('saving to:')
   outfile
   %disp(outfile)
   save(outfile,'an_mat','Sn_mat','t_int','dt','df','inloc',...
                  'Hs_vec','Tp_vec','Tm02_vec');
   %GEN_pause
end

if DO_PLOT
   %%plot time series:
   figure(101);
   subplot(2,1,1);
   plot(tm_vec,Tp_vec,'k');
   hold on;
   plot(tm_vec,Tm02_vec,'--r');
   hold off;
   GEN_proc_fig('time, s','T_p & T_m02, s');
   tstr  = [expt_name,', ' sensor_name];
   ttl   = title(tstr);
   GEN_font(ttl);

   subplot(2,1,2)
   if data_type==2
      Hs_vec   = Hs_vec/2;
      ylab     = 'a_s, ms^{-2}';
   else
      ylab     = 'H_s, m';
   end
   plot(tm_vec,Hs_vec);
   GEN_proc_fig('time, s',ylab);

   %%plot moving window:
   if data_type==2
      ylab1 = 'a, ms^{-2}';
      ylab2 = 'S, m^2s^{-3}';
   else
      ylab1 = 'w, m';
      ylab2 = 'S, m^2s';
   end


   for m=m0:Nwindows
      jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
      time0       = time(jj);
      disp0       = displ(jj);
      disp0       = disp0-mean(disp0);
      %%
      figure(102);
      subplot(1,3,1),plot(time,displ,'k');
      Y  = 1.05*max(displ)*[-1 1];
      ylim(Y);
      hold on;
      plot(time0(1)+0*Y,Y,'r');
      plot(time0(end)+0*Y,Y,'r');
      hold off;
      GEN_proc_fig('time, s',ylab1)
      %%
      subplot(1,3,2);
      plot(time0,disp0,'k');
      ylim(Y);
      xlim([time0(1) time0(end)]);
      GEN_proc_fig('time, s',ylab1)
      ttl   = title(tstr);
      GEN_font(ttl);
      %%
      subplot(1,3,3);
      [Smax,jmax] = max(Sn_mat(:,m));
      periods  = 1./ff;
      Tp = periods(jmax);
      fp = ff(jmax);
      %periods(jmax-1:jmax+1)'

      Sn_max   = max(max(Sn_mat));
      SS       = Sn_max*1.10*[0 1];

      hold off;
      plot(T_target+0*SS,SS,'g')
      hold on;
      plot(periods,Sn_mat(:,m),'k');
      xlim([0 4*T_target]);
      ylim(SS);
      GEN_proc_fig('period, s',ylab2)
      drawnow;
      hold off;
      pause(0.05);
   end
end

function [Sn,ff,an,t_mean,xtra] = get_ft(time,displ,data_type)
%% CALL: [Sn,ff,an,t_mean,xtra] = get_ft(time,displ)
%% time,displ are input vectors
%% data_type is 1 (elevation) or 2 (acceleration - need to integrate F. series wrt time)
%% Sn is the spectrum corresponding to the positive frequencies ff
%%  (ff(1)=0);
%% an are the Fourier coefficients corresponding to [ff(1:end);-fmax;-flipud(ff(end:-1:2))]
%%  (ie no fftshift applied)
%% t_mean   = mean(time)
%% xtra  = {Hs,Tp,T_m02}, T_m02=sqrt(m0/m2);

%% high pass filter by subtracting the mean;
displ = displ-mean(displ);

N        = length(time);
dt       = time(2)-time(1);%%time resolution
T        = time(N)+dt-time(1);%%record length
t_mean   = mean(time);
%%
df    = 1/T;%%freq resolution
fs    = 1/dt;%%sample rate
fmax  = fs/2;%%Nyqvist freq (max possible freq)
%%
ff0   = (-N/2:N/2-1)'*df;
ff    = (0:N/2-1)'*df;
Sn    = 0*ff;
%%
an       = ifft(displ,N);

if 0%data_type==2
   %%was getting strange results from this
   om    = 2*pi*fftshift(ff0);
   an    = an./(-om.^2);
   an(1) = 0;
end

Sn(1)    = abs(an(1)^2)/df;
jj       = 2:N/2;
Sn(jj)   = abs( an(jj).^2 )/df + abs( an(N+1-jj).^2 )/df;
%je       = N/2+1;
%Sn(je)   = 2*abs( an(je).^2 )/df;
%Sn = [abs(an(1))^2;abs(an(2:N2/2)).^2+abs(an(N2:-1:N2/2+2)).^2]/df;

%ff       = (0:df:fmax-df)';
periods  = 1./ff;
%%
[Smax,jmax] = max(Sn);
Tp          = periods(jmax);

m0    = sum(Sn)*df;
m2    = sum(ff.^2.*Sn)*df;
T_m02 = sqrt(m0/m2);%% NB don't need 2\pi factor since moments are wrt freqency.

Hs    = 4*sqrt(m0);
xtra  = {Hs,Tp,T_m02};
%pause

if 0
   tst_var  = [var(displ,1), sum(abs(an.^2)), sum(Sn)*df]
   periods(jmax-1:jmax+1)
%  figure(3)
%  subplot(2,1,1), plot(periods,Sn);
%  %%
%  subplot(2,1,2), plot(time,disp);
%  hold on;
%  disp_app = exp(-2i*pi*time*ff0')*an;
%  plot(time,disp_app,'--r');
%  hold off;
%  pause
end
