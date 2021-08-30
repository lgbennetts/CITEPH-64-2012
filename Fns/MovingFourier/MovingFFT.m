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
% DO_MOV    = create movie 'Mov' if 'vid' is specified in DO_PLOT
% DO_TUNE   = tune the frequency (on by default)
% data_type = 1 (default) wave elevation, 2 acceleration
% data_out
%
% Authors: Timothy Williams & Luke Bennetts
% Begun:   20130722, 09:20:07 CEST

function MovingFFT(file_nm,DO_PLOT,DO_SAVE,DO_DISP)

if ~exist('DO_PLOT','var'); DO_PLOT = 'Aspec-signal'; end
if ~exist('DO_DISP','var'); DO_DISP    = 1;           end
if ~exist('DO_SAVE','var'); DO_SAVE    = 0;           end
if ~exist('DO_MOV','var');  DO_MOV     = 1;           end
if ~exist('DO_TUNE','var'); DO_TUNE    = 1;           end

if ~exist('file_nm','var'); file_nm='eg00001'; end

eval(['load ', file_nm, '-in data description tm Tp data_type '...
 'data_out X Y Tpers wbin'])

if ~exist('data_type','var'); data_type=1; end
if ~exist('Tpers','var'); Tpers = 'Tpers = 20;'; end

eval(Tpers);

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

T_window = Tpers*Tp; %68;
df       = 1/T_window;
Nw       = 2^ceil(log2(T_window/dt));
T_window = dt*Nw;
t_window = (0:dt:(Nw-1))'*dt;

description = [description ';' ' ' num2str(T_window) 's window'];

if strfind(description,'Irregular')
    WaveType ='Irregular';
else
    WaveType ='Regular';
end

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
 if strfind(description,'%')
  dum_i=strfind(description,'%');
  cprintf(0.4*[1,1,1],['>>> ' description(1:dum_i) '%' description(dum_i+1:end) '\n'])
 else
  cprintf(0.4*[1,1,1],['>>> ' description '\n'])
 end
end

%%
Nf_S  = fmax*T_window;
%%
t_fac    = 0.95; % fraction of overlap between windows [1]
t_req    = T_window*(1-t_fac); % time step between windows [seconds]
n_shift  = floor(t_req/dt);
t_req    = n_shift*dt;
t_fac    = 1-t_req/T_window;
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
% Tp_vec   = zeros(Nwindows);
% Hs_vec   = zeros(Nwindows);

for m=m0:Nwindows
 jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
 time0       = tm(jj);
 disp0       = data(jj);
 t_start     = time0(1);
 t_int(m,:)  = time0([1 end]);
 if m==m0 
  [Sn0,tm_vec(m),ff] = get_ft(time0,disp0,data_type);
 else
  [Sn0,tm_vec(m)]    = get_ft(time0,disp0,data_type);
 end
 Sn_mat(:,m) = Sn0; clear Sn0
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots & video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot peak period and sqrt{m0/m2} (av no. down crossings of 0)

fig=fn_getfig;

%%% PEAK PERIOD & SIG WAVE HT %%%
if strfind(DO_PLOT,'TpHm')
 cprintf('green','Tp & Hm not available\n')
 return
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
 
 %  cprintf('red','Paused: hit key to continue...')
 %  pause
 
 if data_type==3
  ylab1 = '\theta, degs';
  ylab2 = 'S, m*degs';
 elseif data_type==2
  ylab1 = 'a, ms^{-2}';
  ylab2 = 'S, m^2s^{-3}';
 else
  ylab1 = 'w [m]';
 end
 
 for m=m0:Nwindows
  jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
  time0       = tm(jj);
  disp0       = data(jj);
  disp0       = disp0-mean(disp0);
  %%
  figure(fig);
  sgtitle(['Fourier Window Video ',description]);
  subplot(1,3,1),plot(tm,data,'k');
  Y0  = 1.05*max(data)*[-1 1];
  ylim(Y0);
  hold on;
  plot(time0(1)+0*Y0,Y0,'r');
  plot(time0(end)+0*Y0,Y0,'r');
  hold off;
  xlabel('t [s]');
  ylabel(ylab1);
%   GEN_proc_fig('t [s]',ylab1)
  %%
  subplot(1,3,2);
  plot(time0,disp0,'k');
  ylim(Y0);
  xlim([time0(1) time0(end)]);
%   GEN_proc_fig('t [s]',ylab1)
  xlabel('t [s]');
  ylabel(ylab1);
  %ttl   = title(tstr);
  %GEN_font(ttl);
  %%
  subplot(1,3,3);
  if findstr('Svid',DO_PLOT)
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
%    GEN_proc_fig('period, s','S, m^2s')
   xlabel('period [s]');
   ylabel('S m^2s');
  elseif findstr('Avid',DO_PLOT)
   [Amax,jmax] = max(an_mat(:,m));
   periods  = 1./ff;
   Tp = periods(jmax);
   fp = ff(jmax);
   an_max   = max(max(an_mat));
   AA       = an_max*1.10*[0 1];
   
   hold off;
   plot(Tp+0*AA,AA,'g')
   hold on;
   plot(periods,an_mat(:,m),'k');
   xlim([0 4*Tp]);
   ylim(AA);
%    GEN_proc_fig('period [s]','amplitude [m]')
   xlabel('period [s]');
   ylabel('amplitude [m]');
  end
  
  drawnow;
  hold off;
  %fn_fullscreen(fig)
  %pause(0.05);
  
  if DO_MOV; Mov(m)=getframe(fig); end
 end
 
 if DO_MOV; save('Movie-MovingFFT','Mov'); end
 
 fig=fig+1;
 
end

%%% ENERGY SPECTRUM %%%
if findstr('Sspec',DO_PLOT)
 figure(fig); fig=fig+1;
 h1=subplot(2,1,1);
 [Smax,jmax] = max(Sn_mat(:,m));
 periods  = 1./ff;
 
 hold off;
%  surf(h1,tm_vec(:,1),periods,log10(Sn_mat))
 surf(h1,tm_vec(:,1),periods,log10(Sn_mat))
 shading interp
 view(h1,[0,90])
 xlim(h1,[tm_vec(1,1), tm_vec(end,1)])
 ylim(h1,[0 4*Tp]);
 ylabel('period, s')
 xlabel('<time>, s')
 %GEN_proc_fig('<time>, s','period, s')
 tstr  = [description '; bin ' int2str(wbin)];
 ttl   = title(tstr);
 %title('log_{10}(S)');
 %GEN_font(ttl)
 cb=colorbar; set(get(cb,'ylabel'),'String','log_{10}(S)'); 
 
 %%% response at harmonics
 
 h2=subplot(2,1,2); hold on
 
 cols = {'b','r','g','c','m','y'};
 
 str = [];
 %Loop 1
%  TpsRegion = [0.8091,0.8296;0.7282,0.7447;0.5958,0.6183];
 TpsRegion = [0.6,1.0];
 %TpsRegion = [0.85,0.95;0.75,0.85;0.65,0.75;0.55,0.65;0.45,0.55];
 for loop=1:size(TpsRegion,1)
%   [~,jj]=min(abs(periods-Tp/loop));
%   plot(h2,tm_vec(:,1),sum(Sn_mat(jj+[-wbin:wbin],:),1),cols{loop})
%   str=[str ' ''\tau_{' int2str(loop) '}=' num2str(periods(jj)) ' '' '  ',' ];
%      plot(h2,tm_vec(:,1),sum(Sn_mat(and(periods <= TpsRegion(loop,1)   , periods >= TpsRegion(loop,2)),:),1),cols{loop})
%     TpI = (TpsRegion(loop,1) + TpsRegion(loop,2))/2;
%     [~,jj]=min(abs(periods-TpI));
%     plot(h2,tm_vec(:,1),sum(Sn_mat(jj+[-wbin:wbin],:),1),cols{loop})

    [~,jjLB]=min(abs(periods-TpsRegion(loop,1)));
    [~,jjUB]=min(abs(periods-TpsRegion(loop,2)));
    LabelName = ['Energy in Period Range (',num2str(round(periods(jjLB),2)),'[s] , ', num2str(round(periods(jjUB),2)) '[s])'];
    plot(h2,tm_vec(:,1),sum(Sn_mat(jjUB:jjLB,:),1),cols{loop},'DisplayName',LabelName);
    str=[str ' ''\tau_{' int2str(loop) '}=' num2str((TpsRegion(loop,1) + TpsRegion(loop,2))/2) ' '' '  ',' ];
 end
 str = str(1:end-1);
 %set(h2,'yscale','log');
 xlabel('<time>, s')
ylabel('S(periods)') 
 %ylabel('log_{10}S(\tau_{i})')
 Y0 = get(h2,'ylim');
 
 legend(h2,'Location','NorthEast');
%  eval(['legend(' str ', ''location'' , ''southeast'' )'])
 
%   h2a = copyobj(h2,gcf);
%  GEN_proc_fig('<time>, s','log_{10}S(\tau_{i})')
 Tind  = fn_Tind(description,Tp,X,tm,WaveType);
 
 TindLB = fn_Tind(description,TpsRegion(1),X,tm,WaveType);
 
 TindUB = fn_Tind(description,TpsRegion(2),X,tm,WaveType);
 if exist('Tind','var')
  hold on
  
  ct_lp=1;
  for loop=1:length(Tind)
   if strfind(Tind(loop).description,'x')
    pp(ct_lp)=plot(h2,TindLB(loop).time + 0*Y0, Y0, [TindLB(loop).colour(1),'--'],'DisplayName',[num2str(TpsRegion(1)) ,'(s)  -',TindLB(loop).description]);
    pp1(ct_lp)=plot(h2,TindUB(loop).time + 0*Y0, Y0, [TindUB(loop).colour(1),':'],'DisplayName',[num2str(TpsRegion(2)) ,'(s)  -',TindUB(loop).description]);

    %ds{ct_lp}=Tind(loop).description; 
    ct_lp=ct_lp+1;
   end
  end
  
  %legend(pp,ds,'Location','NorthEast'); clear pp ds ct_lp
  legend(h2);
  
  hold off
 end % end exist(Tind)
 
 drawnow;
 
end

%%% AMPLITUDE SPECTRUM %%%
if findstr('Aspec',DO_PLOT)
 %an_mat = sqrt(2*Sn_mat);
 
 figure(fig); fig=fig+1;
 h1=subplot(2,1,1);
 [Smax,jmax] = max(an_mat(:,m));
 periods  = 1./ff;
 
 hold off;
 surf(h1,tm_vec(:,1),periods,log10(sqrt(2*Sn_mat)))
 shading interp
 view(h1,[0,90])
 xlim(h1,[tm_vec(1,1), tm_vec(end,1)])
 ylim(h1,[0 4*Tp]);
 xlabel('<time>, s')
 ylabel('period, s')
%  GEN_proc_fig('<time>, s','period, s')
 tstr  = [description '; bin ' int2str(wbin)];
 ttl   = title(tstr);
 %title('log_{10}(S)');
%  GEN_font(ttl)
 cb=colorbar; set(get(cb,'ylabel'),'String','log_{10}(A)'); 
 cb=get(cb,'ylabel'); set(cb,'fontsize',16); clear cb
 
 %%% response at harmonics
 
 h2=subplot(2,1,2); hold on; h2Aspec=h2;
 
 cols = {'b','r','g','c','m','y'};
 
 str = [];
 for loop=1:2
  [~,jj]=min(abs(periods-Tp/loop));
  plot(h2,tm_vec(:,1),sqrt(2*sum(Sn_mat(jj+[-wbin:wbin],:),1)),cols{loop})
  str=[str ' ''\tau_{' int2str(loop) '}=' num2str(periods(jj)) ' '' '  ',' ];
 end
 str = str(1:end-1);
%  set(h2,'yscale','log')
 
 eval(['legend(' str ', ''location'' , ''southeast'' )'])
 xlabel('<time>, s')
 ylabel('A(\tau_{i})')
%  ylabel('log_{10}A(\tau_{i})')
%  GEN_proc_fig('<time>, s','log_{10}A(\tau_{i})')
 
 
 Y0 = get(h2,'ylim');
 
 Tind = fn_Tind(description,Tp,X,tm,WaveType);
 
 if exist('Tind','var')
  hold on
  
  ct_lp=1;
  for loop=1:length(Tind)
   if strfind(Tind(loop).description,'x')
    pp(ct_lp)=plot(h2,Tind(loop).time + 0*Y0, Y0, Tind(loop).colour);
    ds{ct_lp}=Tind(loop).description; ct_lp=ct_lp+1;
   end
  end
  
  %legend(pp,ds,'Location','NorthEast');
  clear pp ds ct_lp
  
 end % end exist(Tind)
 
 drawnow;
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
 
 %clear an_mat
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
 Y0  = 1.05*max(data)*[-1 1];
 ylim(Y0);
% GEN_proc_fig('tm, s',ylab1)
 xlabel(ylab1)
 ylabel(ylab2)
 
 tstr  = description;
 ttl   = title(tstr);
% GEN_font(ttl);
 
 Tind = fn_Tind(description,Tp,X,tm,WaveType);
 
 if exist('Tind','var')
  
  hold on
  ct_lp=1;
  for loop=1:length(Tind)
   if strfind(Tind(loop).description,'x')
    pp(ct_lp)=plot(h1,Tind(loop).time + 0*Y0, Y0, Tind(loop).colour, 'linewidth',2);
    ds{ct_lp}=Tind(loop).description; ct_lp=ct_lp+1;
   end
  end
  
  legend(pp,ds,'Location','NorthEast'); clear pp ds ct_lp
  
 end
 
 drawnow
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%           OUTPUT & SAVE             %%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('data_out','var')
 
 %%% Initialise:
 out_val = []; out_per = []; out_frq = []; out_std = []; out_tim = [];
 out_loc = [];
 
 %%%: Determine time window:
 
 Tind = fn_Tind(description,Tp,X,tm,WaveType);
 tvec=[]; count=1;
 for loop=1:length(Tind)
  if strfind(Tind(loop).description,'x')
   tvec(count)=Tind(loop).time;
   tvecdes{count}=Tind(loop).description; count=count+1;
  end
 end % end loop Tind
 clear count Tind
 [t0,it0]=min(tvec); tvec(it0)=[]; tvecdes(it0)=[]; clear it0
 if strcmp(data_out.tint,'default')
  t1=min(tvec); clear tvec
 else
  eval(data_out.tint)
  clear tvec
 end % end if tint='default'
 tt=find(and(tm_vec>=t0,tm_vec<=t1));
 
 
 %% Loop outputs
 for loop_out=1:length(data_out.name)
  
  %%% Which harmonic:
  if findstr('1st',data_out.name{loop_out})
   harmo=1;
  elseif findstr('2nd',data_out.name{loop_out})
   harmo=2;
  end % end which harmonic
  if ~exist('periods','var'); periods  = 1./ff; end
  [~,jj]=min(abs(periods-Tp/harmo)); pp = periods(jj);
  
  %}{
  %check - mean time, amplitude plots
  % mean amplitude - NSN1 =sum(Sn_mat(jj+[-wbin:wbin],:));
  % plot - plot(tm_vec(:,1),NSN1);
%   maxSnMat = max(max(Sn_mat(:,:)));
%   figure();
%   plot(periods(:),Sn_mat(:,:));
%   hold on;
%   plot([periods(jj- wbin) , periods(jj- wbin) ],[0,1.1*maxSnMat],'--k');
%   plot([periods(jj+ wbin) , periods(jj+ wbin) ],[0,1.1*maxSnMat],'--k');
%   hold off;
%   title('Wave Spectrum At ALL Mean Times')
%   xlabel('Mean Time (s)')
%   ylabel('Wave Spectrum E(t_m)')
% 
%   figure();
%   plot(periods(:),Sn_mat(:,:));
%   hold on;
%   plot([periods(jj- wbin) , periods(jj- wbin) ],[0,1.1*maxSnMat],'--k');
%   plot([periods(jj+ wbin) , periods(jj+ wbin) ],[0,1.1*maxSnMat],'--k');
%   hold off;
%   title('Wave Spectrum At Mean Times In Time Interval')
%   xlabel('Mean Time (s)')
%   ylabel('Wave Spectrum E(t_m)')
%   
%   figure();
%   NSN1 =sum(Sn_mat(jj+[-wbin:wbin],:));
%   maxNSN1 = max(NSN1);
%   plot(tm_vec(:,1),NSN1);
%   hold on;
%   plot([t0,t0],[0,1.1*maxNSN1],'--k');
%   plot([t1,t1],[0,1.1*maxNSN1],'--k');
%   hold off;
% 
%   title('Amplitude at target frequency + neighbours as function of mean times')
%   xlabel('Mean Time (s)')
%   ylabel('Total Amplitude Around Traget Frequency')

  %% Harmonic response
  if findstr('harmo-steady',data_out.name{loop_out})
   
   %%% Amplitude or energy spectra:
   if findstr('amp',data_out.name{loop_out})
    dum_out=sqrt(2*sum(Sn_mat(jj+[-wbin:wbin],tt),1));
   elseif findstr('egy',data_out.name{loop_out})
    dum_out=sum(Sn_mat(jj+[-wbin:wbin],tt));
   end % end amp or ht
  
   %%% OUTPUT %%%
   out_val = [out_val ' ' num2str(mean(dum_out))];
   out_std = [out_std ' ' num2str(std(dum_out))];
   out_per = [out_per ' ' num2str(pp)];
   out_frq = [out_frq ' ' num2str(1/pp)];
   out_tim = [out_tim ' ' '[' num2str(t0) ' ' num2str(t1) ']'];
   out_loc = [out_loc ' ' '[' num2str(X) ' ' num2str(Y) ']'];
   
   %%% PLOT %%%
   if strfind(DO_PLOT,'Aspec')
    plot(h2Aspec,[t0,t1],mean(dum_out)+[0,0],'k','linewidth',5);
   end
   if strfind(DO_PLOT,'signal')
    plot(h1sig,[t0,t1],mean(dum_out)+[0,0],'r','linewidth',5);
   end
   
  %% Time series
  elseif strfind(data_out.name{loop_out},'harmo-time')
   cprintf('green','time series not available\n'); return  
   %%% Amplitude or energy spectra:
   if findstr('amp',data_out.name{loop_out})
    dum_out=0.5*an_mat.*exp(1i*ag_mat);
   elseif findstr('egy',data_out.name{loop_out})
    dum_out=Sn_mat;
   end % end amp or ht
   
   dum_tt{loop_out}     = tm_vec(tt);
   
   %%% TUNE %%%
   if DO_TUNE
    mc   =[tm_vec(tt),1+0*tt]\unwrap(angle(dum_out(jj,tt).'));
    pp   =2*pi*pp/(pp*mc(1) + 2*pi); 
    
    dummer_out{loop_out} = dum_out(jj,tt).*exp(-1i*mc(1)*tm_vec(tt).');
   else
    dummer_out{loop_out} = dum_out(jj,tt);
   end
   
   %%% OUTPUT %%%
   out_val = [out_val ' dummer_out{' int2str(loop_out) '}'];
   out_std = [out_std ' nan'];
   out_per = [out_per ' ' num2str(pp)];
   out_frq = [out_frq ' ' num2str(1/pp)];
   out_tim = [out_tim ' dum_tt{' int2str(loop_out) '}'];
   out_loc = [out_loc ' ' '[' num2str(X) ' ' num2str(Y) ']'];
   
   %%% PLOT %%%
   if strfind(DO_PLOT,'Aspec')
    plot(h2Aspec,tm_vec(tt),2*abs(dum_out(jj,tt)),'k','linewidth',5);
   end
   if strfind(DO_PLOT,'signal')
    plot(h1sig,tm_vec(tt),2*abs(dum_out(jj,tt)),'r','linewidth',5);
   end
   
  end % end output type
  clear dum_out
 end % end loop length(data_out)
 eval(['out_val = {' out_val '};'])
 eval(['out_per = {' out_per '};'])
 eval(['out_frq = {' out_frq '};'])
 eval(['out_std = {' out_std '};'])
 eval(['out_tim = {' out_tim '};'])
 eval(['out_loc = {' out_loc '};'])
 out = struct('value',out_val,'period',out_per,'freq',out_frq,...
  'std',out_std,'time',out_tim,'loc',out_loc);
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
   cprintf('blue',['>>> ' data_out.name{loop_out} ...
    ': ' num2str(out(loop_out).value) ...
    ' (at freq/per ' num2str(out(loop_out).freq) ...
    '/' num2str(out(loop_out).period) ...
    ' std ' num2str(out(loop_out).std) ')\n'])
  end
 end
end % if exist data_out

return

%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%
%%%%%%%%%%%%%%%%%%

% function Tind = fn_Tind(description,Tp,X,tm)

function Tind = fn_Tind(description,Tp,X,tm,WaveType)

if findstr('single floe',description)
  Tind = fn_TestTimes(1/Tp,X,'single',WaveType);
 elseif findstr('calibration',description)
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
 elseif findstr('% conc',description)
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
 elseif strfind(description,'simple harmonic wave')
  dum_t   = {tm(2) tm(end-1)};
  dum_des = {'first waves reach x' 'final waves reach x'};
  dum_col = {'b:' 'r:'};
  Tind = struct('description',dum_des,'time',dum_t,'colour',dum_col);
 end % end if findstr
 
 return

% function [Sn,ff,an,t_mean,xtra] = get_ft(time,displ)
%
% DESCRIPTION: fft in given time window
%              high pass filter by subtracting the mean;
%
% INPUTS:
%
% time      = sampling points in time
% displ     = time series of given quantity
% data_type = 1 (elevation) or 2 (acceleration - need to integrate F. series wrt time)
%             NO LONGER USED!!!!!!
%
% OUTPUTS:
%
% Sn     = spectrum corresponding to the positive frequencies ff (ff(1)=0);
% an     = magnitude of Fourier coefficients corresponding to 
%          [ff(1:end);-fmax;-flipud(ff(end:-1:2))] (ie no fftshift applied)
% ag     = angle corresponding to an
% t_mean = mean(time)
% xtra   = {Hs,Tp,T_m02}, T_m02=sqrt(m0/m2);
% 
% DETAILS:
%
% If displ(n) = eta(tn) for 1 <= n <= N
% &  an(k)    = a(fk)   for 1 <= k <= N,
%
% and we express
%
%            N
% eta(tn) = sum a(fk)*exp(-2i*pi*tn*fk), 1 <= n <= N.
%           k=1
%
%                N
% => displ(n) = sum a(fk)*exp(-2i*pi*(n-1)*(k-1)/N), 1 <= n <= N.
%               k=1
%
% where tn = time(1) + (n-1)*dt is discrete time  
% and   fk =           (k-1)*df is discrete frequency (df=1/N/dt)
%
% then
%
%                N
% a(fk) = (1/N) sum displ(tn)*exp(2i*pi*(n-1)*(k-1)/N), 1 <= k <= N.
%               n=1
%
% L Bennetts & T Williams 2013 / La Seyne & Adelaide

function [Sxx,t_mean,ff] = get_ft(time,displ,data_type)

displ = displ-mean(displ);

N        = length(time);

dt       = time(2)-time(1);                           %%time resolution
T        = time(N)+dt-time(1);                        %%record length
t_mean   = mean(time);
%%
df    = 1/T;                                          %%freq resolution
fs    = 1/dt;                                         %%sample rate
fmax  = fs/2;                                         %%Nyqvist freq (max possible freq)
%%
ff0   = (-N/2:N/2-1)'*df;
ff    = (0:N/2)'*df;
jj    = 1:N/2;
%%
an    = ifft(displ,N);
 
Sxx   = abs(an).^2;

Sneg  = [0;flipud(Sxx(jj+N/2))];                      %%-1,...,-N/2
Spos  = [Sxx(jj);0];                                  %%0,...N/2-1

Sxx   = Spos+Sneg;
 
ff(1) = []; Sxx(1) = [];

%an = sqrt(2*Sxx);

periods  = 1./ff;

%%

if 0
 tst_var  = [var(displ,1), sum(abs(an.^2)), sum(Sn)*df]
 pause
end

return
