function [an,Sn] = citeph_1sensor(time,displ,T_target,outloc,inloc)
%% citeph_1sensor.m
%% Author: Timothy Williams
%% Date:   20130722, 09:20:07 CEST

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
   dir3        = outloc{3};
   outdir   = [outdir0 '/' dir3];
   if ~exist(outdir)
      eval(['!mkdir ' outdir]);
   end
   %%
   sensor_name = outloc{4};%%eg S1,...,S20, A1x, A1y, A1z,...,A6x, A6y, A6z
   outfile     = [outdir,'/singleFFT_',sensor_name,'.mat'];

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
N     = 2^floor(log2(N));
time  = time(end+1-N:end);
displ = displ(end+1-N:end);
[Sn,ff,an,tm,xtra]  = get_ft(time,displ,data_type);

Hs    = xtra{1};
Tp    = xtra{2};
Tm02  = xtra{3};

if DO_SAVE
   df = ff(2)-ff(1);
   dt = time(2)-time(1);
   %disp('saving to:')
   outfile
   %disp(outfile)
   save(outfile,'an','dt','df','inloc',...
                  'Hs','Tp','Tm02');
   %GEN_pause
end

if DO_PLOT

   %%plot moving window:
   if data_type==2
      ylab1 = 'a, ms^{-2}';
      ylab2 = 'S, m^2s^{-3}';
   else
      ylab1 = 'w, m';
      ylab2 = 'S, m^2s';
   end

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
   hold off;
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
