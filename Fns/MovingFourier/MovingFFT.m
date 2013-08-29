% citeph_1sensor_movingFFT.m
% Author: Timothy Williams
% Date:   20130722, 09:20:07 CEST
%
% INPUTS
%
% run_num = specifies runs to analyse via file containing data
%
% INDEP VARIABLES:
%
% tm       = vector
% displ    = vector (same size as tm)
% Tp       = predicted wave period
% outloc   = cell {directory, test, sensor type, sensor} 
% inloc    = string (where the data has come from, to be inc in output)
% 
% DEP VARIABLES:
%
% an_mat 
% Sn_mat
% t_int 
%
% FLAGS:
%
% DO_PLOT
% DO_SIGNAL
% DO_VID
% DO_SAVE
% data_type   = 1 (default) wave elevation, 2 acceleration
% DEL         = delete and saved data after use

function MovingFFT(run_num)

DO_SIGNAL     = 0;
DO_PLOT       = 1;
DO_VID        = 0;
DO_SAVE       = 0;

if ~exist('run_num','var'); run_num=1; end
if ~exist('DEL','var'); DEL=1; end

for run=run_num  

if run > 99
 eval(['load s13',int2str(run),' data description tm Tp data_type X'])
elseif run > 9
 eval(['load s130',int2str(run),' data description tm Tp data_type X'])
else
 eval(['load s1300',int2str(run),' data description tm Tp data_type X'])
end

if ~exist('data_type','var'); data_type=1; end

%cprintf('blue',['>>> ' description '\n'])
disp(['>>> ' description])


% Restrict to one probe
if size(data,2)~=1
 data = data(:,1);
 description = [description ' probe ?!?']; 
end

N     = length(tm);
dt    = tm(2)-tm(1);
T     = tm(end)-tm(1);
fs    = 1/dt;
fmax  = fs/2;
%figure(1),plot(tm,disp);
%%

T_window = 10*Tp;
df       = 1/T_window;
Nw       = 2^ceil(log2(T_window/dt));
T_window = dt*Nw;
t_window = (0:dt:(Nw-1))'*dt;
%%
Nf_S  = fmax*T_window;
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
   time0       = tm(jj);
   disp0       = data(jj);
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
   dt = tm(2)-tm(1);
   %disp('saving to:')
   outfile
   %disp(outfile)
   save(outfile,'an_mat','Sn_mat','t_int','dt','df','inloc',...
                  'Hs_vec','Tp_vec','Tm02_vec');
   %GEN_pause
end

%% Plots & video

%%% plot peak period and sqrt{m0/m2} (av no. down crossings of 0)

fig=101;

if DO_PLOT
   if 0
   %%plot time series:
   figure(fig);fig=fig+1;
   subplot(2,1,1);
   plot(tm_vec,Tp_vec,'k');
   hold on;
   plot(tm_vec,Tm02_vec,'--r');
   hold off;
   GEN_proc_fig('time, s','T_p & T_m02, s');
   %tstr  = [expt_name,', ' sensor_name];
   tstr  = description;
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
   end
end

%%% plot moving window:

if DO_VID 
 
   cprintf('red','Paused: hit key to continue...')
   pause
 
   if data_type==2
      ylab1 = 'a, ms^{-2}';
      ylab2 = 'S, m^2s^{-3}';
   else
      ylab1 = 'w, m';
      ylab2 = 'S, m^2s';
   end

   for m=m0:Nwindows
      jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
      time0       = tm(jj);
      disp0       = data(jj);
      disp0       = disp0-mean(disp0);
      %%
      figure(fig); fig=fig+1;
      subplot(1,3,1),plot(tm,data,'k');
      Y  = 1.05*max(data)*[-1 1];
      ylim(Y);
      hold on;
      plot(time0(1)+0*Y,Y,'r');
      plot(time0(end)+0*Y,Y,'r');
      hold off;
      GEN_proc_fig('tm, s',ylab1)
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
      plot(Tp+0*SS,SS,'g')
      hold on;
      plot(periods,Sn_mat(:,m),'k');
      xlim([0 4*Tp]);
      ylim(SS);
      GEN_proc_fig('period, s',ylab2)
      drawnow;
      hold off;
      pause(0.05);
   end
   
%%% spectogram   
   
elseif DO_PLOT 
 
   figure(fig); fig=fig+1;
   h1=subplot(2,1,1);
   [Smax,jmax] = max(Sn_mat(:,m));
   periods  = 1./ff;
   
   hold off;
   surf(h1,tm_vec(:,1),periods,log10(Sn_mat))
   shading interp
   view(h1,[0,90])
   xlim(h1,[tm_vec(1,1), tm_vec(end,1)]) 
   ylim(h1,[0 4*Tp]);
   GEN_proc_fig('time, s','period, s')
   tstr  = description;
   ttl   = title(tstr);
   %title('log_{10}(S)');
   GEN_font(ttl) 
   colorbar
   
%%% response at harmonics   
   
   h2=subplot(2,1,2); hold on  
   
   cols = {'b','r','g','c','m','y'};
   
   str = [];
   for loop=1:2
    [~,jj]=min(abs(periods-Tp/loop));
    plot(h2,tm_vec(:,1),Sn_mat(jj,:),cols{loop})
    str=[str ' ''\tau_{' int2str(loop) '}=' num2str(periods(jj)) ' '' '  ',' ];
   end
   str = str(1:end-1);
   set(h2,'yscale','log')
   
   eval(['legend(' str ', ''location'' , ''southeast'' )'])
   
   GEN_proc_fig('time, s','log_{10}S(\tau_{i})')
   
   Y = get(h2,'ylim');
   
   Tind = fn_TestTimes(1/Tp,X);
   
   hold on
   
   for loop=6:length(Tind)
    plot(h2,Tind(loop).time + 0*Y, Y,'k:');
   end
   
   hold off
   
   drawnow;
   
end % end DO_VID/DO_PLOT

%%% signal

if and(DO_SIGNAL,~DO_VID)
 if data_type==2
  ylab1 = 'a, ms^{-2}';
  ylab2 = 'S, m^2s^{-3}';
 else
  ylab1 = 'w, m';
  ylab2 = 'S, m^2s';
 end
   
   figure(fig); fig=fig+1;
   plot(tm,data,'k');
   Y  = 1.05*max(data)*[-1 1];
   ylim(Y);
   GEN_proc_fig('tm, s',ylab1)
   
   tstr  = description;
   ttl   = title(tstr);
   GEN_font(ttl);

   hold on
   
   Tind = fn_TestTimes(1/Tp,X);
   
   for loop=6:length(Tind)
    plot(Tind(loop).time + 0*Y, Y,'r--');
   end
   
   %hold off
end

end % end run

return

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

function [Sn,ff,an,t_mean,xtra] = get_ft(time,displ,data_type)

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

if 0
   Sn(1)    = abs(an(1)^2)/df;
   jj       = 2:N/2;
   Sn(jj)   = abs( an(jj).^2 )/df + abs( an(N+1-jj).^2 )/df;
else
   jj = 1:N/2;
   aneg  = flipud(an(jj+N/2));%%-1,...,-N/2
   apos  = an(jj);%%0,...N/2-1
   Sn    = (abs(apos.^2)+abs(aneg.^2))/df;
   ff    = ff+df/2;%%move output frequency to middle of bin instead of the left;
end
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
   pause
end

return
