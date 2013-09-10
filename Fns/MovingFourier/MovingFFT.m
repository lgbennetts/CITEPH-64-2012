% function MovingFFT(file_nm,DO_PLOT,DO_SAVE)
%
% INPUTS
%
% file_pre = string that locates the file
% run_num  = specifies runs to analyse via file containing data
%
% VARIABLES:
%
% tm       = raw time
% tm_vec   = mid-pt of Fourier windows
% N        = number of points in series
% dt       = time step
% T        = time span
% fs       = sampling frequency
% fmax     = fs/2
% T_pers   = appx number of periods to fit in window
% T_window = size of Fourier window
% df       = 1/T_window
% Nw       = number of points in time window as power of 2
% t_window = times in window shifted to begin at 0
% Nf_S     = fmax*T_window
% n_shift  = floor(t_req/dt);
% t_req    = n_shift*dt
% Nwindows = number of windows
% Sn_mat   = energy spectrum as function of period and tm_vec
% an_mat   = amplitude spectrum as function of period and tm_vec
%
% FLAGS:
%
% DO_PLOT   = 'Spec'   wave energy spectrum
%             'Aspec'  amplitude spectrum
%             'vid'    video of energy spectrum
%             'signal' raw signal
%             'TpHm'   peak period & sig wave height spectra
% DO_DISP   = display output (on by default)
% DO_SAVE   = save data (off by default)
% data_type = 1 (default) wave elevation, 2 acceleration
% data_out
%
% Authors: Timothy Williams & Luke Bennetts
% Begun:   20130722, 09:20:07 CEST


function MovingFFT(file_nm,DO_PLOT,DO_SAVE,DO_DISP)

if ~exist('DO_PLOT','var'); DO_PLOT = 'Aspec-signal'; end
if ~exist('DO_DISP','var'); DO_DISP    = 1;           end
if ~exist('DO_SAVE','var'); DO_SAVE    = 0;           end

if ~exist('file_nm','var'); file_nm='s00001'; end

eval(['load ', file_nm, '-in data description tm Tp data_type '...
 'data_out X T_pers'])

if ~exist('data_type','var'); data_type=1; end
if ~exist('T_pers','var'); T_pers = 20; end

% Restrict to one probe
if size(data,2)~=1
 cprintf('green','>>> Restricting to 1 data set\n');
 data = data(:,1);
end

N     = length(tm);
dt    = tm(2)-tm(1);
T     = tm(end)-tm(1);
fs    = 1/dt;
fmax  = fs/2;

%%

T_window = T_pers*Tp;
df       = 1/T_window;
Nw       = 2^ceil(log2(T_window/dt));
T_window = dt*Nw;
t_window = (0:dt:(Nw-1))'*dt;

description = [description ';' ' ' num2str(T_window) 's window'];
  
if DO_SAVE
 if and(exist('data_out','var'),~exist([file_nm '.mat']))
  descriptions{1} = description;
  eval(['save ' file_nm ' descriptions;'])
  clear descriptions
  sv_num=1;
 elseif exist('data_out','var')
  eval(['load ' file_nm ' descriptions;'])
  sv_num=length(descriptions)+1;
  descriptions{sv_num}=description;
  eval(['save ' file_nm ' descriptions -append;'])
  clear descriptions
 end
end

if DO_DISP
 cprintf(0.4*[1,1,1],['>>> ' description '\n'])
end

%%
Nf_S  = fmax*T_window;
%%
t_req    = 0.5;%s
n_shift  = floor(t_req/dt);
t_req    = n_shift*dt;
%%
Nwindows = floor((T-T_window)/t_req);
Sn_mat   = zeros(Nf_S,Nwindows);
an_mat   = zeros(Nf_S,Nwindows);

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

%% SAVE 1

% if and(XTRA_SV,exist('data_out','var'))
% %  df = ff(2)-ff(1);
% %  dt = tm(2)-tm(1);
% %  save(file_nm,'an_mat','Sn_mat','t_int','dt','df',...
% %   'Hs_vec','Tp_vec','Tm02_vec','T_pers','-append');
% end

clear T_pers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots & video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot peak period and sqrt{m0/m2} (av no. down crossings of 0)

fig=101;

%%% PEAK PERIOD & SIG WAVE HT %%%
if strfind(DO_PLOT,'TpHm')
 if 0
  %%plot time series:
  figure(fig);fig=fig+1;
  subplot(2,1,1);
  plot(tm_vec,Tp_vec,'k');
  hold on;
  plot(tm_vec,Tm02_vec,'--r');
  hold off;
  GEN_proc_fig('<time>, s','T_p & T_m02, s');
  %tstr  = [expt_name,', ' sensor_name];
  tstr  = description;
  ttl   = title(tstr);
  GEN_font(ttl);
  
  subplot(2,1,2)
  if data_type==3
   ylab     = '\theta_s, degs';
  elseif data_type==2
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

if findstr('vid',DO_PLOT)
 
 cprintf('red','Paused: hit key to continue...')
 pause
 
 if data_type==3
  ylab1 = '\theta, degs';
  ylab2 = 'S, m*degs';
 elseif data_type==2
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
  GEN_proc_fig('<time>, s',ylab1)
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
 
 %%% ENERGY SPECTRUM %%%
elseif findstr('Sspec',DO_PLOT)
 
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
 GEN_proc_fig('<time>, s','period, s')
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
 
 GEN_proc_fig('<time>, s','log_{10}S(\tau_{i})')
 
 Y = get(h2,'ylim');
 
 if findstr('single floe',description)
  Tind = fn_TestTimes(1/Tp,X,'single');
 elseif findstr('% conc',description)
  Tind = fn_TestTimes(1/Tp,X,'attn');
 end % end if findstr
 
 if exist('Tind','var')
  hold on
  
  for loop=1:length(Tind)
   if strfind(Tind(loop).description,'x')
    plot(h2,Tind(loop).time + 0*Y, Y, Tind(loop).colour);
   end
  end
  
  hold off
 end % end exist(Tind)
 
 drawnow;
 
 %%% AMPLITUDE SPECTRUM %%%
elseif findstr('Aspec',DO_PLOT)
 
 figure(fig); fig=fig+1;
 h1=subplot(2,1,1);
 [Smax,jmax] = max(an_mat(:,m));
 periods  = 1./ff;
 
 hold off;
 surf(h1,tm_vec(:,1),periods,log10(an_mat))
 shading interp
 view(h1,[0,90])
 xlim(h1,[tm_vec(1,1), tm_vec(end,1)])
 ylim(h1,[0 4*Tp]);
 GEN_proc_fig('<time>, s','period, s')
 tstr  = description;
 ttl   = title(tstr);
 %title('log_{10}(S)');
 GEN_font(ttl)
 colorbar
 
 %%% response at harmonics
 
 h2=subplot(2,1,2); hold on; h2Aspec=h2;
 
 cols = {'b','r','g','c','m','y'};
 
 str = [];
 for loop=1:2
  [~,jj]=min(abs(periods-Tp/loop));
  plot(h2,tm_vec(:,1),an_mat(jj,:),cols{loop})
  str=[str ' ''\tau_{' int2str(loop) '}=' num2str(periods(jj)) ' '' '  ',' ];
 end
 str = str(1:end-1);
 set(h2,'yscale','log')
 
 eval(['legend(' str ', ''location'' , ''southeast'' )'])
 
 GEN_proc_fig('<time>, s','log_{10}A(\tau_{i})')
 
 Y = get(h2,'ylim');
 
 if findstr('single floe',description)
  Tind = fn_TestTimes(1/Tp,X,'single');
 elseif findstr('% conc',description)
  Tind = fn_TestTimes(1/Tp,X,'attn');
 end % end if findstr
 
 if exist('Tind','var')
  hold on
  
  for loop=1:length(Tind)
   if strfind(Tind(loop).description,'x')
    plot(h2,Tind(loop).time + 0*Y, Y, Tind(loop).colour, 'linewidth',2);
   end
  end

 end % end exist(Tind)
 
 drawnow;
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
 
end % end DO_VID/DO_PLOT

%%% RAW SIGNAL %%%

if findstr('signal',DO_PLOT)
 if data_type==3
  ylab1 = '\theta, degs';
  ylab2 = 'S, m*degs';
 elseif data_type==2
  ylab1 = 'a, ms^{-2}';
  ylab2 = 'S, m^2s^{-3}';
 else
  ylab1 = 'w, m';
  ylab2 = 'S, m^2s';
 end
 
 figure(fig); fig=fig+1;
 h1 = subplot(1,1,1); h1sig = h1;
 plot(h1,tm,data,'k');
 Y  = 1.05*max(data)*[-1 1];
 ylim(Y);
 GEN_proc_fig('tm, s',ylab1)
 
 tstr  = description;
 ttl   = title(tstr);
 GEN_font(ttl);
 
 if findstr('single floe',description)
  Tind = fn_TestTimes(1/Tp,X,'single');
 elseif findstr('% conc',description)
  Tind = fn_TestTimes(1/Tp,X,'attn');
 end % end if findstr
 
 if exist('Tind','var')
  
  hold on
  
  for loop=1:length(Tind)
   if strfind(Tind(loop).description,'x')
    plot(h1,Tind(loop).time + 0*Y, Y, Tind(loop).colour, 'linewidth',2);
   end
  end
  
 end
 
 drawnow
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
 
end

%% OUTPUT & SAVE

if exist('data_out','var')
 out_val = []; out_per = []; out_frq = []; out_std = [];
 for loop_out=1:length(data_out)
  if findstr('harmo-steady',data_out(loop_out).name)
   if findstr('amp',data_out(loop_out).name)
    dum_out=an_mat;
   elseif findstr('egy',data_out(loop_out).name)
    dum_out=Sn_mat;
   end % end amp or ht
   if findstr('1st',data_out(loop_out).name)
    harmo=1;
   elseif findstr('2nd',data_out(loop_out).name)
    harmo=2;
   end % end which harmonic
   if findstr('single floe',description)
    Tind = fn_TestTimes(1/Tp,X,'single');
   elseif findstr('% conc',description)
    Tind = fn_TestTimes(1/Tp,X,'attn');
   end % end get test times
   if ~exist('periods','var'); periods  = 1./ff; end
   [~,jj]=min(abs(periods-Tp/harmo)); pp = periods(jj);
   t_vec=[]; count=1;
   for loop=1:length(Tind)
    if strfind(Tind(loop).description,'x')
     t_vec(count)=Tind(loop).time; count=count+1;
    end
   end % end loop Tind
   clear count Tind
   [t0,it0]=min(t_vec); t_vec(it0)=[]; clear it0
   if strcmp(data_out(loop_out).tint,'default')
    t1=min(t_vec); clear t_vec
   else
    eval(data_out(loop_out).tint)
    clear t_vec
   end % end if tint='default'
   tt=find(and(tm_vec>=t0,tm_vec<=t1));
   %%% OUTPUT %%%
   out_val = [out_val ' ' num2str(mean(dum_out(jj,tt)))];
   out_std = [out_std ' ' num2str(std(dum_out(jj,tt)))];
   out_per = [out_per ' ' num2str(pp)];
   out_frq = [out_frq ' ' num2str(1/pp)];
   %%% PLOT %%%
   if strfind(DO_PLOT,'Aspec')
    plot(h2Aspec,[t0,t1],mean(dum_out(jj,tt))+[0,0],'k--','linewidth',5);
   end 
   if strfind(DO_PLOT,'signal')
    plot(h1sig,[t0,t1],mean(dum_out(jj,tt))+[0,0],'r:','linewidth',5);
   end
  end % end output type
  clear dum_out
 end % end loop length(data_out)
 eval(['out_val = {' out_val '};'])
 eval(['out_per = {' out_per '};'])
 eval(['out_frq = {' out_frq '};'])
 eval(['out_std = {' out_std '};'])
 out = struct('value',out_val,'period',out_per,'freq',out_frq,'std',out_std);
 if DO_SAVE
  if sv_num==1
   outputs{1}=out;
   eval(['save ' file_nm ' outputs -append']);
  else
   eval(['load ' file_nm ' outputs']);
   outputs{sv_num}=out;
   eval(['save ' file_nm ' outputs -append']); clear outputs
  end
 end
 if DO_DISP
  for loop=1:length(out)
   cprintf('blue',['>>> ' data_out(loop_out).name ...
    ': ' num2str(out(loop_out).value) ...
    ' (at freq/per ' num2str(out(loop_out).freq) ...
    '/' num2str(out(loop_out).period) ...
    ' std ' num2str(out(loop_out).std) ')\n'])
  end
 end
end % if exist data_out

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
 %%%%%%%%% TW (July '13) %%%%%%%%%
 %  ff    = (0:N/2-1)'*df;
 %  Sn    = 0*ff;
 %  aneg  = flipud(an(jj+N/2));%%-1,...,-N/2
 %  apos  = an(jj);%%0,...N/2-1
 %  Sn    = (abs(apos.^2)+abs(aneg.^2))/df;
 %  ff    = ff+df/2;%%move output frequency to middle of bin instead of the left;
 %%%%%%%%% LB (Aug  '13) %%%%%%%%%
 ff    = (0:N/2)'*df;
 Sn    = 0*ff;
 aneg  = [0;flipud(an(jj+N/2))];                      %%-1,...,-N/2
 apos  = [an(jj);0];                                  %%0,...N/2-1
 Sn    = (abs(apos.^2)+abs(aneg.^2))/df;
 ff(1) = []; Sn(1) = [];
 apos(1) = []; aneg(1) = [];
 an    = 2*abs(apos);       % nb. aneg = conj(apos) for a real signal
 %%%%%%%%% END LB VS TW  %%%%%%%%%
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
