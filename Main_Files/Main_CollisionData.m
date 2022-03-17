% function Main_CollisionData
%
% DESCRIPTION:
%
% INPUTS:
%
% DO_CALC       = if need to do (or redo) calculation set to 1;
%                 else (if it's done) set to 0;
% DO_ANALYSIS   = if it's done (& the file's saved) can speed things up by
%                 setting this to 0;
% DO_PLOT       = plot data (1) or not (0);
% DO_DISP       = display info (1) or not (0);
% DO_SAVE       = save data (1) or not (0);
% conc          = concentration (set to 79 to get collisions)
% opts          = which accelerometers, e.g. opts = {'Ax','Ay','Az'};
% xtra_opts     = e.g. 'front', 'middle', 'back' (rows)
% data_out      = specify analysis window 
% fig           = figure handle
% col           = string for color, linestyle, markerstyle
%
% accelerometers: front row  (A3,A6)
%                 middle row (A4,A5)
%                 back row   (A1,A2)
%
% T Williams Oct 2013 / Bergen
%
% REVISION HISTORY:
%                  modified by L Bennetts Oct 2013 / Adelaide

function Main_CollisionData

%% PRELIMINARIES:

if ~or(strcmp(getenv('LOGNAME'),'a1612881'),...
  strcmp(getenv('LOGNAME'),'lbennetts'))
 OTHER_USR=1;
end

if exist('OTHER_USR','var')
 DO_CALC = 1;
 DO_ANALYSIS = 1;
 DO_SAVE = 1;
 conc=79;
end

if ~exist('fig','var');  fig=fn_getfig; end
if ~exist('col','var');  col=' ''b.'' , ''markersize'' , 12'; end
if ~exist('conc','var'); conc=79; end
if ~exist('HT','var');   HT=fn_WhatTestData(conc,'Irregular',0); 
                         HT=HT(1:2,:); end
if ~exist('DO_CALC','var');     DO_CALC     = 1; end
if ~exist('DO_ANALYSIS','var'); DO_ANALYSIS = 1; end
if ~exist('DO_PLOT','var');     DO_PLOT     = 1; end
if ~exist('DO_DISP','var');     DO_DISP     = 1; end
if ~exist('DO_SAVE','var');     DO_SAVE     = 1; end
if ~exist('opts','var');        opts        = {'Ax'}; end
if ~exist('xtra_opts','var');   xtra_opts   = 'middle'; end
if ~exist('data_out','var');
 % data_out.tint='ts=ts+4*T_target; tf=ts+10*T_target; ';
 % data_out.tint='tf=min(t_vec); ';
 data_out.tint=["Tpers=fn_Tpers(T_target,'Irregular');[ts,tf] = fn_tint(t_vec,ts,T_target,Tpers);"];
end

if conc==39
 test_type   = 'c39';
elseif conc==79
 test_type   = 'c79';
end

%% CALCULATIONS, ETC:

citeph_coll0(test_type,HT,opts,xtra_opts,data_out, ...
 DO_CALC,DO_ANALYSIS,DO_PLOT,DO_DISP,DO_SAVE,fig,col)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function citeph_coll0(test_type,HT,opts,xtra_opts,data_out,...
 DO_CALC,DO_ANALYSIS,DO_PLOT,DO_DISP,DO_SAVE,fig,col)

tres  = .004;     %%time resolution [s]

%% Find tests:

%% c79
if strcmp(test_type,'c79')
 %  Ntests      = 19;
 %  Ntests_reg  = 15;
 if isempty(HT)
  %   Ntests      = 12;
  %   Ntests_reg  = 8;
  %   test_list   = 1:Ntests_reg;
  HT=fn_WhatTestData(79,'Irregular',0); HT=HT(1:2,end);
  test_list=1:size(HT,2);
 else
  eval(['c_prams = conc79_testspecs();'])
  pers = zeros(1,length(c_prams)); hts = pers;
  for loop=1:length(c_prams)
   if strcmp(c_prams(loop).type,'Irregular')
    pers(loop)=10\c_prams(loop).period;
    hts(loop) =10*c_prams(loop).wave_height;
   end
  end
  test_list=[];
  for lp=1:size(HT,2)
   dum_test = find(and(pers==HT(2,lp),hts==HT(1,lp)));
   if length(dum_test)>1
    cprintf('magenta',['>>> test repeated ' int2str(length(dum_test)) ...
     ' times\n'])
    if 0
     pers(dum_test(2:end))=[];
     hts(dum_test(2:end)) =[];
     dum_test=dum_test(1);
     cprintf('magenta',['... using 1st test only\n'])
    else
     cprintf('magenta',['... using all ' int2str(length(dum_test)) ...
      ' tests\n'])
    end
   end
   test_list=[test_list, dum_test]; clear dum_test
  end
 end
 %% c39
elseif strcmp(test_type,'c39')
 if isempty(HT)
  %   Ntests      = 12;
  %   Ntests_reg  = 8;
  %   test_list   = 1:Ntests_reg;
  HT=fn_WhatTestData(39,'Regular',0); HT=HT(1:2,:);
  test_list=1:size(HT,2);
 else
  eval(['c_prams = conc39_testspecs();'])
  pers = zeros(1,length(c_prams)); hts = pers;
  for loop=1:length(c_prams)
   if strcmp(c_prams(loop).type,'Regular')
    pers(loop)=10\c_prams(loop).period;
    hts(loop) =10*c_prams(loop).wave_height;
   end
  end
  test_list=[];
  for lp=1:size(HT,2)
   dum_test = find(and(pers==HT(2,lp),hts==HT(1,lp)));
   if length(dum_test)>1
    cprintf('magenta',['>>> test repeated ' int2str(length(dum_test)) ...
     ' times\n'])
    dum_test=dum_test(1);
    cprintf('magenta',['... using 1st test only\n'])
   end
   test_list=[test_list, dum_test]; clear dum_test
  end
  clear pers hts
 end
elseif strcmp(test_type,'calib')
 Ntests   = 17;
end

%% SENSOR LOOP

for n=1:length(opts)
 opt      = opts{n};
 if strcmp(test_type,'c79')
  outfile   = ['Collision_data/Conc79/coll_stats_' opt '-' xtra_opts '.mat'];
  figfile   = ['Collision_data/Conc79/collision_strength_' opt '-' xtra_opts '.fig'];
  figfileb  = ['Collision_data/Conc79/collision_hist_' opt '-' xtra_opts '.fig'];
 elseif strcmp(test_type,'c39')
  outfile   = ['Collision_data/Conc39/coll_stats_' opt '-' xtra_opts '.mat'];
  figfile   = ['Collision_data/Conc39/collision_strength_' opt '-' xtra_opts '.fig'];
  figfileb  = ['Collision_data/Conc39/collision_hist_' opt '-' xtra_opts '.fig'];
 end
 
 %% CALCULATIONS: Get accelerometer data: 
 
 if DO_CALC
  if DO_DISP; cprintf(0.4*[1,1,1],'Beginning calculations...\n'); end
  ct_test=1;
  for test_num=test_list
   cprintf(0.4*[1,1,1],['>>> experiment: H=' num2str(hts(test_num)) ...
    '; T=' num2str(pers(test_num)) '\n']); ct_test=ct_test+1;
   citeph_coll_1test(test_type,test_num,opt,xtra_opts,data_out,...
    DO_PLOT,DO_DISP,or(DO_ANALYSIS,DO_SAVE))
  end
  clear ct_test
  if DO_DISP; cprintf(0.4*[1,1,1],'... ending calculations\n'); end
  
 end
 
 %% ANALYSIS: Extract quantities of interest:
 
 if DO_ANALYSIS==1
  if DO_DISP; cprintf(0.4*[1,1,1],'Beginning analysis...\n'); end
  ct_test = 1;
  T_coll=zeros(length(test_list),1); col_dt = T_coll;
  rms_a = T_coll; isat = rms_a; Hi_Frac = rms_a; c_freq = rms_a;
  for test_num=test_list
   %%open results & plot;
   [Tcol,Acol,Asig,Tp,Hs,sensor_names,Col_dt,isat0,Cfq] =...
    citeph_coll_1test_OpenResultsFile(test_type,test_num,opt,xtra_opts,...
    DO_CALC);
   T_target(ct_test,1) = Tp;
   H_target(ct_test,1) = Hs;
   if DO_DISP; 
    cprintf(0.4*[1,1,1],['>>> T=' num2str(Tp) '; H=' num2str(Hs) '\n'])
   end
   Nsens = length(Tcol); r_ct=0;
   for r=1:Nsens 
    if ~isempty(Tcol{r})
     r_ct = r_ct+1;
     acol2 = Acol{r};
     %%collision period (mean length of time between collisions);
     tcol2               = Tcol{r};
     dt                  = tcol2(2:end)-tcol2(1:end-1);
     T_coll(ct_test)   = T_coll(ct_test)+mean(dt);%%6x1 vector
     dt_all{ct_test}   = dt;
     %%collision freq per period
     c_freq(ct_test)   = c_freq(ct_test) + Cfq{r};
     %%
     col_dt(ct_test)   = col_dt(ct_test) + mean(Col_dt{r});
     %%
     rms_a(ct_test)    = rms_a(ct_test) + mean(abs(acol2)); %sqrt(mean(acol2.^2));%%6x1 vector
     isat(ct_test)     = isat(ct_test)+max(isat0{r});
     %%fraction of variance that is high frequency (percentage)
     asig1             = Asig(r,1);
     asig2             = Asig(r,2);
     Hi_Frac(ct_test)  = Hi_Frac(ct_test) + (asig2/asig1)^2*100;%%6x1 vector
    end
   end
   T_coll(ct_test)  = T_coll(ct_test)/r_ct; 
   col_dt(ct_test)  = col_dt(ct_test)/r_ct;
   c_freq(ct_test)  = c_freq(ct_test)/r_ct;
   rms_a(ct_test)   = rms_a(ct_test)/r_ct; 
   Hi_Frac(ct_test) = Hi_Frac(ct_test)/r_ct;
   isat(ct_test)    = isat(ct_test)~=0; clear r_ct
   if ~exist('out')
    eval(['!mkdir out']);
   end
   %plot(T_target,rms/H_target,'x');
   %hold on
   ct_test = ct_test + 1;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SAVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if DO_CALC
   if DO_SAVE
    eval(['!cp -r Collision_data/temp/ Collision_data;']);
    eval(['!rm -r Collision_data/temp/;']);
   else
    eval(['!rm -r Collision_data/temp/;']);
   end
  end
  if DO_SAVE
   save(outfile,'rms_a','T_target','H_target','sensor_names',...
     'T_coll','Hi_Frac','col_dt','isat','xtra_opts'); %,'dt_all'
  end 
  if DO_DISP; cprintf(0.4*[1,1,1],'... ending analysis\n'); end
 elseif DO_ANALYSIS==-1%%do/don't do calc
  load(outfile);
 end
end%n - opt (Ax,Ay,Az)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if and(DO_PLOT,abs(DO_ANALYSIS))
 figure(fig);
 Nsens = size(rms_a,2);
 %%  COLLISION STRENGTHS:
 subplot(2,1,1); hold on; set(gca,'box','on','fontsize',16)
 title('Root mean squared acceleration vs period')
 for n=1:Nsens
  eval(['P(n)  = plot(T_target,rms_a(:,n),' col ');']) %./H_target
  if ~isempty(find(isat(:,n)))
   Q(n)  = plot(T_target(find(isat(:,n))),rms_a(find(isat(:,n)),n),'rs'); %...
   %./H_target(find(isat(:,n))),'rs');
  end
  [~,IA] = unique(T_target); IA = [0;IA];
  jj_nl=[]; dum_vec=find(diff(IA)>1);
  cnt_ct=1;
  for ct=1:length(dum_vec)
   loop=dum_vec(ct);
   dum_inds=IA(loop)+1:IA(loop+1);
   if length(unique(H_target(dum_inds)))>1
    [~,jj_nl(cnt_ct)]=max(H_target(dum_inds));
    jj_nl(cnt_ct)=dum_inds(jj_nl(cnt_ct)); cnt_ct=cnt_ct+1;
   end
  end
  clear ct IA dum_inds loop cnt_ct
  if ~isempty(jj_nl)
   R(n)  = plot(T_target(jj_nl),rms_a(jj_nl,n),'r.','markersize',12); %...
   %./H_target(jj_nl),'r.','markersize',12);
  end
 end
 if and(exist('Q','var'),exist('R','var'))
  legend([P,R,Q],{['A-' xtra_opts],'large amp','saturation'},'location','northwest');
 elseif exist('Q','var')
  legend([P,Q],{['A-' xtra_opts],'saturation'},'location','northwest');
 elseif exist('R','var')
  legend([P,R],{['A-' xtra_opts],'large amp'},'location','northwest');
 else
  legend(P,{['A-' xtra_opts]},'location','northwest');
 end
 hold off;
 %GEN_proc_fig('T_p, s','rms(a_{col})/H_s, s^{-2}');
 xlabel('T_p, s','fontsize',12)
 ylabel('<a_{col}>, ms^{-2}','fontsize',12) %/H_s
 %% COLLISION FREQUENCY:
 subplot(2,1,2); hold on; set(gca,'box','on','fontsize',16)
 title('Collision frequency per period vs period','fontsize',16)
 for n=1:Nsens
  eval(['P(n)  = plot(T_target,c_freq(:,n),' col ');'])
  if ~isempty(jj_nl)
   R(n)  = plot(T_target(jj_nl),c_freq(jj_nl,n)...
    ,'r.','markersize',12);
  end
 end
 %legend(P,sensor_names,'location','eastoutside');
 %ylim([0 3.5]);
 hold off;
 %GEN_proc_fig('T_p, s','T_{col}/T_p');
 xlabel('T_p, s','fontsize',12)
 ylabel('c_{fq}/T_p','fontsize',12)
%  subplot(2,1,2); hold on;
%  for n=1:Nsens
%   eval(['P(n)  = plot(T_target,T_coll(:,n)./T_target,' col ');'])
%   if ~isempty(jj_nl)
%    R(n)  = plot(T_target(jj_nl),T_coll(jj_nl,n)./T_target(jj_nl)...
%     ,'r.','markersize',12);
%   end
%  end
%  %legend(P,sensor_names,'location','eastoutside');
%  %ylim([0 3.5]);
%  hold off;
%  %GEN_proc_fig('T_p, s','T_{col}/T_p');
%  xlabel('T_p, s','fontsize',12)
%  ylabel('T_{col}/T_p','fontsize',12)
 %%
%  subplot(3,1,3); hold on;
%  for n=1:Nsens
%   P(n)  = plot(T_target,Hi_Frac(:,n),col); %colin{n});
%  end
%  %legend(P,sensor_names,'location','eastoutside');
%  hold off;
%  %GEN_proc_fig('T_p, s','High-freq variance %');
%  xlabel('T_p, s','fontsize',12)
%  ylabel('High-freq variance %','fontsize',12)
 if DO_SAVE; saveas(gcf,figfile,'fig'); end
 fn_fullscreen
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HISTOGRAMS:

if and(0*DO_PLOT,abs(DO_ANALYSIS))
 figure(fig);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 loc   = [5 6 1 3 4 2];
 for r=1:Nsens
  %for r=4
  N1          = size(dt_all,1);
  Dt_all      = [];
  check_ind   = 1;
  for test_num=1:N1
   Tp       = T_target(test_num);
   Hs       = H_target(test_num);
   %Dt_all   = [Dt_all;dt_all{test_num,r}/Tp];
   if check_ind
    Dt_all   = dt_all{test_num,r};
    Y        = Dt_all/Tp;
    Yav      = mean(Y);
    figure(fig)
    subplot(3,2,loc(r));
    hist(Y,50);
    ttl   = title([sensor_names{r},': ',num2str(Tp),' s, ',num2str(Hs),'m; mean = ',num2str(Yav,'%0.1f')]);
    GEN_font(ttl);
    hold on;
    yl = get(gca,'ylim');
    plot(Yav+0*yl,yl,'g');
    ylim(yl);
    hold off;
    GEN_proc_fig('\Delta t/T_p','frequency');
   else%% bunch them together
    Dt_all   = [Dt_all;dt_all{test_num,r}];
   end
  end
  % {r,loc(r),sensor_names{r}}
  if ~check_ind
   figure(fig);
   subplot(3,2,loc(r));
   hist(Dt_all,50);
   ttl   = title(sensor_names{r});
   GEN_font(ttl);
   GEN_proc_fig('\Delta t/T_p','frequency');
  end
  drawnow
  fn_fullscreen;
 end
 if DO_SAVE; saveas(gcf,figfileb,'epsc'); end
end

return

% function citeph_coll_1test(test_type,test_num,opt)

function citeph_coll_1test(test_type,test_num,opt,xtra_opts,data_out,...
 DO_PLOT,DO_DISP,DO_SAVE)

if ~exist('cls_fig','var'); cls_fig=1; end

if nargin==0
 test_type   = 'c79';
 test_num    = 19;
end

[time,data,file_list,sensor_names]  = citeph_get_data(test_type,test_num,opt,'true');
% nsens = length(file_list);

inds = fn_inds(test_type,opt,xtra_opts);

time=time(:,inds); data=data(:,inds);
file_list=file_list(inds); sensor_names=sensor_names(inds);

if 0
 figure(1);data(:,nsens)
 plot(time(:,6),data(:,nsens))
 pause
end

if strcmp(test_type,'c79')
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num,'true');
 %tt2   = 'conc_79/';
 tt2 = 'Conc79/';
else%% 'c39'
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num,'true');
 %tt2   = 'conc_39/';
 tt2 = 'Conc39/';
end

Nsens = size(data,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basedir  = citeph_user_specifics;
% if basedir(end)=='/'; basedir(end)=[]; end
% jj       = find(basedir=='/');
% basedir2 = [basedir,'/../PROCESSED_data/'];
if DO_SAVE
 basedir2 = 'Collision_data/';
 if ~exist(basedir2)
  eval(['mkdir ' basedir2]);
 end
 basedir2 = [basedir2 'temp/'];
 if ~exist(basedir2)
  eval(['mkdir ' basedir2]);
 end
 basedir2 = [basedir2,tt2];
 if ~exist(basedir2)
  eval(['mkdir ' basedir2]);
 end
else
 basedir2=[];
end
% if strcmp(type,'Regular')
%  type  = 'regular';
% else
%  type  = 'irregular';
% end
% basedir3 = [basedir2,type,'/'];%%/regular or /irregular
% if ~exist(basedir3)
%  eval(['mkdir ' basedir3]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figs=[];

if DO_DISP; cprintf(0.4*[1,1,1],['>>> expt num:   ' expt_dir '\n']); end

for n=1:Nsens
 %for n=nsens
 time0 = time(:,n);
 displ = data(:,n);
 inloc = file_list{n};
 if DO_DISP; cprintf(0.4*[1,1,1],['>>> sensor: ' sensor_names{n} '\n']); end
 
 %%
 outloc                  = {basedir2,expt_dir,opt,sensor_names{n} };
 %[an_mat,Sn_mat,t_int]   = citeph_1sensor_collisions(time0,displ,T_target,outloc,inloc);
 dumfig = citeph_1sensor_collisions(time0,displ,T_target,H_target,outloc,...
  inloc,data_out,DO_PLOT,DO_DISP,DO_SAVE);
 
 figs=[figs, dumfig]; clear dumfig
end

if cls_fig; close(figs); end

% if DO_PLOT
%  pause
%  for lp=1:length(figs); close(figure(figs(lp))); end
% end

return

% citeph_coll_1test_OpenResultsFile

function [Tcol,Acol,Asig,T_target,H_target,sensor_names,Col_dur,isat,cfq] =...
 citeph_coll_1test_OpenResultsFile(test_type,test_num,opt,xtra_opts,DO_CALC)

if nargin==0
 test_type   = 'c79';
 test_num    = 19;
end

[time,data,file_list,sensor_names]  = citeph_get_data(test_type,test_num,opt,'true');
% nsens = 6;

inds = fn_inds(test_type,opt,xtra_opts);

time=time(:,inds); data=data(:,inds);
file_list=file_list(inds); sensor_names=sensor_names(inds);

if 0
 figure(1);data(:,nsens)
 plot(time(:,6),data(:,nsens))
 pause
end

if strcmp(test_type,'c79')
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num,'true');
 %tt2   = 'conc_79/';
 tt2 = 'Conc79/';
else%% 'c39'
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num,'true');
 %tt2   = 'conc_39/';
 tt2 = 'Conc39/';
end

Nsens = size(data,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basedir  = citeph_user_specifics;
% if basedir(end)=='/'; basedir(end)=[]; end
% jj       = find(basedir=='/');
% basedir2 = [basedir,'/../PROCESSED_data/'];
basedir2='Collision_data/';
if DO_CALC; basedir2=[basedir2 'temp/']; end
expt_name   = ['T' sprintf('%3.2f',T_target) '_H' sprintf('%3.2f',H_target)];
basedir3 = [basedir2,tt2,expt_name,'_',expt_dir,'/',opt,'/mat_files/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:Nsens
 time0 = time(:,n);
 displ = data(:,n);
 inloc = file_list{n};
 fname = [basedir3,'collisions_',sensor_names{n},'.mat'];
 if exist(fname,'file')
  load(fname);
  Tcol{n}     = tcol2;
  Acol{n}     = acol2;
  Asig(n,:)   = [asig1,asig2];
  Col_dur{n}  = col_dt;
  isat{n}     = isat2;
  cfq{n}      = c_freq;
 else
  cprintf('magenta',['>>> ' fname ' does not exist\n']) 
  Tcol{n}     = [];
  Acol{n}     = [];
  Asig(n,:)   = [0,0];
  Col_dur{n}  = [];
  isat{n}     = [];
 end
end

return

%% citeph_get_data.m
%% Author: Timothy Williams
%% Date:   20130723, 12:21:29 CEST
%% CALL: [time,data] = citeph_get_data(test_type,test_num)
%% test_type:  'c79', 'c39', 'calib'
%% test_num:   1,2,...
%% opt:        'S', 'S_zoom', 'Ax', 'Ay', 'Az'
%% scale:      'true' (measurement/basin scale), 'full' (lengths increased by 100, times increased by 10)

function [time,data,file_list,sensor_names] = citeph_get_data(test_type,test_num,opt,scale)

if ~exist('what_tests','var'); what_tests = 'final'; end

if nargin==0
 %test_type   = 'calib';
 %test_num    = 13;
 %opt         = 'S';
 test_type   = 'c79';
 test_num    = 1;
 %  opt         = 'S';
 %  opt         = 'S_zoom';
 %  opt         = 'Ax';
 %  opt         = 'Ay';
 opt         = 'Az';
end
if ~exist('scale')
 scale = 'false';%%full scale is default
end
basedir  = citeph_user_specifics(what_tests);

if basedir(end)=='/'
 basedir(end)=[];
end

if strcmp(test_type,'calib');%%real or calibration
 if strcmp(what_tests,'prelim')
  fdir  = [basedir '/calibration/data_Wave_Calibration/calibration_waves/'];
 elseif strcmp(what_tests,'final')
  fdir  = [basedir '/1-Calibration/calibration_accelerometers/'];
 end
 str1  = 'calib_houle_';
 %%
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_calib_prams(test_num);
 %%
 if strcmp(type,'Regular');%% reg/irreg waves
  str2  = 'reg_';
 else
  str2  = 'irreg_';
 end
elseif strcmp(test_type,'c79')
 if strcmp(what_tests,'prelim')
  fdir  = [basedir '/attenuation_tests/Conc79/'];
 elseif strcmp(what_tests,'final')
  fdir  = [basedir '/2-Wave_attenuation/Conc79/Irregular/'];
 end
 str1  = 'houle_';
 %%
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num);
 %%
 if strcmp(type,'Regular');%% reg/irreg waves
  str2  = 'reg_';
 else
  str2  = 'irr_';
 end
else%% 'c39'
 if strcmp(what_tests,'prelim')
  fdir  = [basedir '/attenuation_tests/Conc39/'];
 elseif strcmp(what_tests,'final')
  fdir  = [basedir '/2-Wave_attenuation/Conc39/'];
 end
 str1  = 'houle_';
 %%
 [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num);
 %%
 if strcmp(type,'Regular');%% reg/irreg waves
  str2  = 'reg_';
 else
  str2  = 'irr_';
 end
end

if strcmp(test_type,'calib');
 %% prev 20 = full eta series
 %% last 20 = zoom eta series
 nn    = [20 20];
 if strcmp(opt,'S')
  sub_dir  = '/wave_probes/';
 else
  sub_dir  = '/wave_probes_zoom/';
 end
 nstep = 1;
else
 %% prev 20 = full eta series
 %% prev 20 = zoom eta series
 %% last 18 = [a_x,a_y,a_z] series
 if strcmp(test_type,'c79')
  nn = [20 20 16];
 elseif strcmp(test_type,'c39')
  nn = [20 20 14];
 end
 if strcmp(opt,'S')
  %sub_dir  = '/wave_probes/';
  n0 = 1:20;
  for j=1:20
   sensor_names{j}   = ['S',num2str(j)];
  end
 elseif strcmp(opt,'S_zoom')
  %sub_dir  = '/wave_probes_zoom/';
  n0 = 21:40;
  nstep = 1;
  for j=1:20
   sensor_names{j}   = ['S',num2str(j)];
  end
 elseif strcmp(opt,'Ax')
  %sub_dir  = '/accelerometers/';
  if strcmp(test_type,'c79')
   n0 = 41:56;
   nvec  = [2 4 6 8 11 14];
   jj=1:6;
   for j=jj
    sensor_names{j}   = ['A',num2str(j),'x'];
   end
  elseif strcmp(test_type,'c39')
   n0 = 41:54;
   nvec  = [2 4 6 9 12];
   jj=[1,2,4:6]; j_ct=1;
   for j=jj
    sensor_names{j_ct}   = ['A',num2str(j),'x']; j_ct=j_ct+1;
   end
   clear j_ct
  end
  
 elseif strcmp(opt,'Ay')
  %sub_dir  = '/accelerometers/';
  if strcmp(test_type,'c79')
   n0 = 41:56;
   nvec  = [1 5 7 9 12 15];
   jj = 1:6;
   for j=jj
    sensor_names{j}   = ['A',num2str(j),'y'];
   end
  elseif strcmp(test_type,'c39')
   n0 = 41:54;
   nvec  = [1 5 7 10 13];
   jj=[1,2,4:6]; j_ct=1;
   for j=jj
    sensor_names{j_ct}   = ['A',num2str(j),'x']; j_ct=j_ct+1;
   end
   clear j_ct
  end
 else%% 'Az'
  %Nfiles,nn(3)
  %sub_dir  = '/accelerometers/';
  if strcmp(test_type,'c79')
   n0 = 41:56;
   nvec  = [3 10 13 16];
   jz    = [1,4:6];%%sensors that had the z-acceleration;
   jj = 1:4;
   for j=jj
    sensor_names{j}   = ['A',num2str(jz(j)),'z'];
   end
  elseif strcmp(test_type,'c39')
   n0 = 41:54;
   nvec  = [3 8 11 14];
   jj=[1,2,4:6]; j_ct=1;
   for j=jj
    sensor_names{j_ct}   = ['A',num2str(j),'x']; j_ct=j_ct+1;
   end
   clear j_ct
  end
 end
end

%%

fx_dir   = [fdir,expt_dir]; %,sub_dir];
DD       = dir([fx_dir ,'/', str1,'*.dat']);
n0       = length(DD)-sum(nn)+n0;
DD       = DD(n0); clear n0
Nprobes  = length(DD);

if strcmp(opt,'S')%%get full time series from all probes
 file_list   = cell(Nprobes,1);
 for j=1:Nprobes
  fname          = [fx_dir '/' DD(j).name];
  file_list{j}   = fname;
  %%
  A  = load(fname);
  if j==1
   time  = A(:,1);
   N     = length(time);
   data  = zeros(N,Nprobes);
  end
  {data,A}
  data(:,j)   = A(:,2);
 end
elseif strcmp(opt,'S_zoom')%%zoomed timeseries
 file_list   = cell(Nprobes,1);
 for j=1:Nprobes
  fname          = [fx_dir '/' DD(j).name];
  file_list{j}   = fname;
  A              = load(fname);
  %%
  if j==1
   N     = size(A,1);
   time  = zeros(N,Nprobes);
   data  = zeros(N,Nprobes);
  end
  %%
  time(:,j)   = A(:,1);
  data(:,j)   = A(:,2);
 end
else%%accelerometer full
 Nprobes     = length(nvec);
 file_list   = cell(Nprobes,1);
 for j=1:Nprobes
  fname          = [fx_dir '/' DD(nvec(j)).name];
  file_list{j}   = fname;
  A              = load(fname);
  %%
  if j==1
   N     = size(A,1);
   time  = zeros(N,1);
   data  = zeros(N,Nprobes);
  end
  %%
  time(:,j)   = A(:,1);
  data(:,j)   = A(:,2);
 end
end

if strcmp(scale,'true')
 fac      = 100;
 time     = time/sqrt(fac);
 %data     = data/fac;
end

return

%% citeph_get_calib_prams.m
%% Author: Timothy Williams
%% Date:   20130724, 10:14:53 CEST

function [dirname,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num,scale)

DO_TEST  = 0;
if nargin==0
 test_num = 19;
 DO_TEST  = 1;
elseif 0
 test_num = 13;
 DO_TEST  = 1;
end

if ~exist('scale')
 scale = 'false';
end
c79_prams   = conc79_testspecs();

dirname     = c79_prams(test_num).dirname;
expt_name   = c79_prams(test_num).name;
T_target    = c79_prams(test_num).period;
H_target    = c79_prams(test_num).wave_height;
type	    = c79_prams(test_num).type;

% if strcmp(type,'Regular')
%    dirname  = ['regular/' dirname];
% else
%    dirname  = ['irregular/' dirname];
% end

if strcmp(scale,'true')%%scale variables to measured/basin scale;
 T_target = T_target/10;
 H_target = H_target/100;
end

if DO_TEST
 dirname,T_target,H_target,type,expt_name
end

return

%% citeph_get_calib_prams.m
%% Author: Timothy Williams
%% Date:   20130724, 10:14:53 CEST

function [dirname,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num,scale)

DO_TEST  = 0;
if nargin==0
 test_num = 13;
 DO_TEST  = 1;
end

if ~exist('scale')
 scale = 'false';
end
c39_prams = conc39_testspecs();

dirname     = c39_prams(test_num).dirname;
expt_name   = c39_prams(test_num).name;
T_target    = c39_prams(test_num).period;
H_target    = c39_prams(test_num).wave_height;
type	    = c39_prams(test_num).type;

% if strcmp(type,'Regular')
%    dirname  = ['regular/' dirname];
% else
%    dirname  = ['irregular/' dirname];
% end

if strcmp(scale,'true')%%scale variables to measured/basin scale;
 T_target = T_target/10;
 H_target = H_target/100;
end

if DO_TEST
 dirname,T_target,H_target,type,expt_name
end

return

function inds = fn_inds(test_type,opt,xtra_opts)

inds = [];

 %%%%%%%%%%%%%%%
 %% c79 %%%%%%%%
 %%%%%%%%%%%%%%%
if strcmp(test_type,'c79')
 %%% ALL %%%
 if strfind(xtra_opts,'all')
  if strfind(xtra_opts,'left')
   if strfind(opt,'z')
    inds=2;   %% A4
   else
    inds=[3,4,2];
   end
  elseif strfind(xtra_opts,'right')
   if strfind(opt,'z')
    inds=[4,3,1];   %% A6,A5,A1
   else
    inds=[6,5,1];
   end
  else
   if strfind(opt,'z')
    inds=[4,3,2,1];     %% A1,A4,A5,A6
   else
    inds=[6,3,5,4,1,2];
   end
  end
 end
 %%% FRONT %%%
 if strfind(xtra_opts,'front')
  if strfind(xtra_opts,'front-left')
   if strfind(opt,'z')
    cprintf('magenta',['A3z does not exist\n'])
   else
    inds=[inds,3];
   end
  elseif strfind(xtra_opts,'front-right')
   if strfind(opt,'z')
    inds=[inds,4];    %% A6
   else
    inds=[inds,6];
   end
  else
   if strfind(opt,'z')
    inds=[inds,4];
   else
    inds=[inds,3,6];
   end
  end
 end
 %%% BACK %%%
 if strfind(xtra_opts,'back')
  if strfind(xtra_opts,'back-left')
   if strfind(opt,'z')
    cprintf('magenta',['A2z does not exist\n'])
   else
    inds=[inds,2];
   end
  elseif strfind(xtra_opts,'back-right')
   inds=[inds,1];
  else
   inds=[inds,1,2];
  end
 end
 %%% MIDDLE %%%
 if strfind(xtra_opts,'middle')
  if strfind(xtra_opts,'middle-left')
   if strfind(opt,'z')
    inds=[inds,4];
   else
    inds=[inds,4];
   end
  elseif strfind(xtra_opts,'middle-right')
   if strfind(opt,'z')
    inds=[inds,4];
   else
    inds=[inds,5];
   end
  else
   if strfind(opt,'z')
    inds=[inds,2,3];
   else
    inds=[inds,4,5];
   end
  end
 end
 %%%%%%%%%%%%%%%
 %% c39 %%%%%%%%
 %%%%%%%%%%%%%%%
elseif strcmp(test_type,'c39')
 %%% ALL %%%
 if strfind(xtra_opts,'all')
  if strfind(xtra_opts,'left')
   if strfind(opt,'z')
    inds=2;
   else
    inds=[2,3];
   end
  elseif strfind(xtra_opts,'right')
   if strfind(opt,'z')
    inds=[4,3,1];
   else
    inds=[5,4,1];
   end
  else
   if strfind(opt,'z')
    inds=[4,3,2,1];
   else
    inds=[5,4,3,2,1];
   end
  end
 end
 %%% FRONT %%%
 if strfind(xtra_opts,'front')
  if strfind(xtra_opts,'front-left')
    cprintf('magenta',['A5 does not exist\n'])
  elseif strfind(xtra_opts,'front-right')
   if strfind(opt,'z')
    inds=4;
   else
    inds=5;
   end
  else
   if strfind(opt,'z')
    inds=[inds,4];
   else
    inds=[inds,5];
   end
  end
 end
 %%% BACK %%%
 if strfind(xtra_opts,'back')
  if strfind(xtra_opts,'back-left')
   if strfind(opt,'z')
    cprintf('magenta',['A2z does not exist\n'])
   else
    inds=2;
   end
  elseif strfind(xtra_opts,'back-right')
   inds=1;
  else
   if strfind(opt,'z')
    inds=[inds,1];
   else
    inds=[inds,1,2];
   end
  end
 end
 %%% MIDDLE %%%
 if strfind(xtra_opts,'middle')
  if strfind(xtra_opts,'middle-left')
   if strfind(opt,'z')
    inds=[inds,2];
   else
    inds=[inds,3];
   end
  elseif strfind(xtra_opts,'middle-right')
   if strfind(opt,'z')
    inds=[inds,3];
   else
    inds=[inds,4];
   end
  else
   if strfind(opt,'z')
    inds=[inds,2,3];
   else
    inds=[inds,3,4];
   end
  end
 end
end

return

%% citeph_1sensor_movingFFT.m
%% Author: Timothy Williams
%% Date:   20130722, 09:20:07 CEST

function fig=citeph_1sensor_collisions(time,displ,T_target,...
 H_target,outloc,inloc,data_out,DO_PLOT,DO_DISP,DO_SAVE)

if max(abs(displ))==0
 cprintf('magenta',['Warning -> no data for T=' num2str(T_target) ...
  ' H=' num2str(H_target) ' ' outloc{4} '\n'])
 fig=[];
 return
end

% Saturation of signal:
if ~exist('sat','var'); sat=20; end

if ~exist('DO_SAVE','var'); DO_SAVE     = 0; end
data_type   = 1;%%default (1) is wave elevation, 2->acceleration
%%
expt_name      = '';
sensor_name    = '';

if DO_PLOT; fig_ct=0; else; fig=[]; end

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

sensor_name = outloc{4};%%eg S1,...,S20, A1x, A1y, A1z,...,A6x, A6y, A6z
 
if DO_PLOT
 if strcmp(sensor_name(3),'x')
  ylab  = 'a_x, ms^{-2}';
 elseif strcmp(sensor_name(3),'y')
  ylab  = 'a_y, ms^{-2}';
 elseif strcmp(sensor_name(3),'z')
  ylab  = 'a_z, ms^{-2}';
 end
end

if and(DO_SAVE,exist('outloc'))
 %DO_SAVE  = 1;
 %%
 dir1  = outloc{1};
 if ~exist(dir1)
  eval(['!mkdir ' dir1]);
 end
 %%
 expt_name   = ['T' sprintf('%3.2f',T_target) '_H' sprintf('%3.2f',H_target)];
 expt_name   = [expt_name '_' outloc{2}];%%date,time,year;
 %outdir0     = dir1;
 outdir0     = [dir1 '/' expt_name];
 if ~exist(outdir0)
  eval(['!mkdir ' outdir0]);
 end
 %%
 dir3     = outloc{3};%% eg Az
 outdir   = [outdir0 '/' dir3];
 if ~exist(outdir)
  eval(['!mkdir ' outdir]);
 end
 %%
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
 figfile4    = [outdir2,'/collisions_Aspec_',sensor_name,'.eps'];
 %%
 %figfile2    = [outdir,'/collisions_',sensor_name,'.mat'];
 
 if strcmp(sensor_name(1),'A')
  data_type   = 2;
 end
end

N     = length(time);
%%
dt    = time(2)-time(1);
T     = time(end)-time(1);
fs    = 1/dt;               % sampling frequency
fmax  = fs/2;
%figure(1),plot(time,disp);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
 tstart      = 50;%%wave maker starts [s]
 dist        = 13.244;%%dist to wave maker [m]
 om          = 2*pi/T_target;
 waveno      = om^2/9.81;
 groupvel    = om/waveno/2;%%group velocity [m/s] is half the phase velocity at infinite depth
 trav_time   = dist/groupvel;
 ts          = tstart+trav_time;
 %if DO_DISP; cprintf(0.4*[1,1,1],['>>> start time: ' num2str(ts) 's\n']); end
 tstop = T_target*120;%%approx time wave maker stops
 tf    = tstop+trav_time;
else
 if strfind(outloc{4},'3')
  X = -2;
 elseif strfind(outloc{4},'6')
  X = -2;
 elseif strfind(outloc{4},'4')
  X = 0;
 elseif strfind(outloc{4},'5')
  X = 0;
 elseif strfind(outloc{4},'1')
  X = 2;
 elseif strfind(outloc{4},'2')
  X = 2;
 end
 Tind = fn_TestTimes(1/T_target,X,'attn','Irregular'); clear X
 count=1;
 for loop=1:length(Tind)
  if strfind(Tind(loop).description,'x')
   t_vec(count)=Tind(loop).time; count=count+1;
  end
 end % end loop Tind
 clear count Tind
 [ts,it0]=min(t_vec); t_vec(it0)=[]; clear it0
 if strcmp(data_out.tint,'default')
  tstop=min(t_vec); clear t_vec
 else
  eval(data_out.tint)
  ts = min(t_vec);
  tf = max(t_vec);
  clear t_vec
 end % end if tint='default'
end

j_rel = find(time>=ts & time<=tf);%%when wave-maker starts
if mod(length(j_rel),2)==1
 j_rel(1) = [];
end
t_rel    = time(j_rel);
a_rel0   = displ(j_rel);
a_rel    = a_rel0-mean(a_rel0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FREQ SPECTRUM:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_filt   = T_target/4;%%s
f_filt   = 1/T_filt;%%Hz
[Bh,Ah]  = butter(3,f_filt/fmax,'high');
[Bl,Al]  = butter(3,f_filt/fmax,'low');

%% Filter, FFT & Detect Collisions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter

ah = filter(Bh,Ah,a_rel);
al = filter(Bl,Al,a_rel);

%%% FFT

[Sn_a ,ff,an,t_mean,xtra]    = citeph_get_ft(t_rel,a_rel); %T_target,xtra
asig1 = xtra{1}/2;
Tp1   = xtra{2}; clear xtra
[Sn_ah,~,an_h,~,xtra]        = citeph_get_ft(t_rel,ah);
asig2    = xtra{1}/2;
Tp2   = xtra{2}; clear xtra
[Sn_al,~,an_l,~,xtra]        = citeph_get_ft(t_rel,al);
asig3 = xtra{1}/2;
Tp3   = xtra{2}; clear xtra

%%% Detect collisions

a_thresh = 2*asig3;                    %% 4 std deviations in amp of acc  %0.0287; %
jcol     = find(abs(ah)>a_thresh);     %% collision indices
tcol     = t_rel(jcol);                %% collision times
acol     = ah(jcol);                   %% collision amplitudes
isat     = abs(acol)>=sat;            %% saturation reached
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time between big accelerations is too small => same collision
if 1
 Tc    = 10*dt;%%max collision time
 %%
 tcol2a         = [];
 tcol2b         = [];
 tcol_current   = [];
 acol2          = [];
 isat2          = [];
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
   if ~isempty(find(acol(jc)>=sat))        %% Saturation reached
    isat2=[isat2;1];
   else
    isat2=[isat2;0];
   end
   acol2 = [acol2; mean(abs(acol(jc)))]; %sqrt( mean(acol(jc).^2) )];
   %%
   tcol_current   = [];
  end
  %tcol_current,pause
 end
 
 tcol2    = (tcol2a+tcol2b)/2;
 col_dt   = tcol2a-tcol2b;
else
 tcol2    = tcol;
 acol2    = acol;
 col_dt   = dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency of collisions per period

c_freq = length(acol2)*T_target/(t_rel(end)-t_rel(1));

%% Save
if DO_SAVE==1
 save(outfile,'tcol2','acol2','col_dt','asig1','asig2','asig3',...
  'T_target','H_target','isat2','sat','c_freq');
end

clear c_freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Raw time series + high/low pass filters

if and(DO_PLOT,1)
 %yl = [-10 21];
 xl=[t_rel(1),t_rel(end)];
 fig_ct=fig_ct+1;
 fig(fig_ct)=fn_getfig();
 figure(fig(fig_ct));
 subplot(3,1,1), plot(t_rel,a_rel);
 j0       = strfind(inloc,'Conc');
 inloc0   = inloc(j0:end);
 inloc0   = strrep(inloc0,'_','\_');
 
 ttl   = {['Raw time series (',sensor_name, '): a_s = ',num2str(asig1,'%0.2e'),' ms^{-2}; ',...
  'T_p = ',num2str(T_target), ' s; H_s = ',num2str(H_target), ' m'],inloc0};
 %GEN_font(ttl);
 title(ttl,'fontsize',12)
 xlim(xl); %ylim(yl);
 %GEN_proc_fig('Time, s',ylab)
 xlabel('Time, s','fontsize',12)
 ylabel(ylab,'fontsize',12)
 %%
 subplot(3,1,3), plot(t_rel,al);
 ttl   = ['Low-pass-filtered time series: ' ...
  ,num2str(100*(asig3/asig1)^2,'%0.1f'),' % of variance'];
 %GEN_font(ttl);
 %GEN_proc_fig('Time, s',ylab)
 title(ttl,'fontsize',12)
 xlabel('Time, s','fontsize',12)
 ylabel(ylab,'fontsize',12)
 xlim(xl);
 %%
 subplot(3,1,2);
 plot(t_rel,ah);%%plot once 1st to get y-range
 %%
 yl = get(gca,'ylim'); %a_thresh*[-1,1];
 for n=1:length(tcol2)
  if n==1
   hold on;
  end
  plot(tcol2(n)+0*xl,yl,'g:');
  if isat2(n); plot(tcol2(n),log10(sat),'ro','markersize',12); end
 end
 hold on;
 plot(xl,0*xl+a_thresh,'m');
 plot(xl,0*xl-a_thresh,'m');
 %%plot(t_rel,ah);%%plot again so it's on top of other lines
 hold off;
 ttl   = ['High-pass-filtered time series: ' ...
  ,num2str(100*(asig2/asig1)^2,'%0.1f'),' % of variance' ...
  '; threshold ' sprintf('%3.2f',a_thresh)];
 xlim(xl); %ylim(yl);
 %GEN_font(ttl);
 title(ttl,'fontsize',12)
 %GEN_proc_fig('Time, s',ylab)
 xlabel('Time, s','fontsize',12)
 ylabel(ylab,'fontsize',12)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 if DO_SAVE
  saveas(gcf,figfile1,'epsc');
 end
 
 fn_fullscreen;
 
end

%% plot histogram:

if and(DO_PLOT,0)
 fig_ct=fig_ct+1;
 fig(fig_ct)=fn_getfig();
 figure(fig(fig_ct));
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
 %set(gca,'yscale','log','xscale','log');
 GEN_proc_fig(xlab2,'Cumulative probability');
 
 if nr==2
  subplot(nr,2,3);
  plot(X2,Y2);
  GEN_proc_fig(xlab2b,'Probability density');
  
  subplot(nr,2,4);
  yy = dx2*fliplr(cumsum(fliplr(Y2)));
  plot(X2,yy);
  %set(gca,'yscale','log','xscale','log');
  GEN_proc_fig(xlab2b,'Cumulative probability');
 end
 
 if DO_SAVE; saveas(gcf,figfile2,'epsc'); end
 
 fn_fullscreen;
end

%% Plot Sn spectra:

if and(DO_PLOT,0)
 fig_ct=fig_ct+1;
 fig(fig_ct)=fn_getfig();
 figure(fig(fig_ct));
 subplot(3,1,1), plot(1./ff,Sn_a);
 set(gca,'yscale','log');
 hold on;
 yl = get(gca,'ylim');
 p(1)=plot(T_target+0*yl,yl,'r');
 p(2)=plot(Tp1+0*yl,yl,'--k');
 p(3)=plot(T_filt+0*yl,yl,'g');
 %xlim([0 15]);
 %GEN_proc_fig('Period, s','SD, m^2s^{-3}');
 xlabel('Period, s','fontsize',12)
 ylabel('SD, m^2s^{-3}','fontsize',12)
 ttl   = title(['Spectral density (',sensor_name,')']);
 legend(p,{'target','peak period','filter'},'Location','SouthEast'); clear p
 %GEN_font(ttl);
 title(ttl,'fontsize',12)
 hold off;
 %%
 subplot(3,1,2), plot(1./ff,Sn_ah);
 [Smax,jmax]  = max(Sn_ah);
 set(gca,'yscale','log');
 hold on;
 yl = get(gca,'ylim');
 p(1)=plot(T_target+0*yl,yl,'r');
 p(2)=plot(Tp2+0*yl,yl,'--k');
 p(3)=plot(T_filt+0*yl,yl,'g');
 xlim([0 2*T_target]);
 %xl = get(gca,'xlim');
 plot(xl,0*xl+Smax,'k')
 hold off;
 %GEN_proc_fig('Period, s','SD (HPF), m^2s^{-3}');
 xlabel('Period, s','fontsize',12)
 ylabel('SD (HPF), m^2s^{-3}','fontsize',12)
 legend(p,{'target','peak period','filter'},'Location','SouthEast'); clear p
 %%
 subplot(3,1,3), plot(1./ff,Sn_al);
 set(gca,'yscale','log');
 hold on;
 yl = get(gca,'ylim');
 p(1)=plot(T_target+0*yl,yl,'r');
 p(2)=plot(Tp3+0*yl,yl,'--k');
 p(3)=plot(T_filt+0*yl,yl,'g');
 %xlim([0 15]);
 hold off;
 %GEN_proc_fig('Period, s','SD (LPF), m^2s^{-3}');
 xlabel('Period, s','fontsize',12)
 ylabel('SD (LPF), m^2s^{-3}','fontsize',12)
 legend(p,{'target','peak period','filter'},'Location','SouthEast'); clear p
 if DO_SAVE
  saveas(gcf,figfile3,'epsc');
 end
 fn_fullscreen;
end

%% Plot an spectra:

if and(DO_PLOT,0)
 fig_ct=fig_ct+1;
 fig(fig_ct)=fn_getfig();
 figure(fig(fig_ct));
 subplot(3,1,1), plot(1./ff,an);
 set(gca,'yscale','log');
 hold on;
 ymn=min(an); ymx=max(an);
 yl = 10.^[floor(log10(ymn)),ceil(log10(ymx))]; clear ymn ymx
 set(gca,'ylim',yl)
 p(1)=plot(T_target+0*yl,yl,'r');
 p(2)=plot(T_filt+0*yl,yl,'g');
 %xlim([0 15]);
 %GEN_proc_fig('Period, s','SD, m^2s^{-3}');
 xlabel('Period, s','fontsize',12)
 ylabel('A, m^2s^{-3}','fontsize',12)
 ttl   = ['Spectral density (',sensor_name,')'];
 legend(p,{'target','filter'},'Location','SouthEast'); clear p
 %GEN_font(ttl);
 title(ttl,'fontsize',12)
 hold off;
 %%
 subplot(3,1,2), plot(1./ff,an_h);
 [Smax,jmax]  = max(Sn_ah);
 set(gca,'yscale','log');
 hold on;
 ymn=min(an_h); ymx=max(an_h);
 yl = 10.^[floor(log10(ymn)),ceil(log10(ymx))]; clear ymn ymx
 set(gca,'ylim',yl)
 p(1)=plot(T_target+0*yl,yl,'r');
 p(2)=plot(T_filt+0*yl,yl,'g');
 %xlim([0 2*T_target]);
 xl = get(gca,'xlim');
 plot(xl,0*xl+Smax,'k')
 hold off;
 %GEN_proc_fig('Period, s','SD (HPF), m^2s^{-3}');
 xlabel('Period, s','fontsize',12)
 ylabel('A (HPF), m^2s^{-3}','fontsize',12)
 legend(p,{'target','filter'},'Location','SouthEast'); clear p
 %%
 subplot(3,1,3), plot(1./ff,an_l);
 set(gca,'yscale','log');
 hold on;
 ymn=min(an_l); ymx=max(an_l);
 yl = 10.^[floor(log10(ymn)),ceil(log10(ymx))]; clear ymn ymx
 set(gca,'ylim',yl)
 p(1)=plot(T_target+0*yl,yl,'r');
 p(2)=plot(T_filt+0*yl,yl,'g');
 %xlim([0 15]);
 hold off;
 %GEN_proc_fig('Period, s','SD (LPF), m^2s^{-3}');
 xlabel('Period, s','fontsize',12)
 ylabel('A (LPF), m^2s^{-3}','fontsize',12)
 legend(p,{'target','filter'},'Location','SouthEast'); clear p
 if DO_SAVE
  saveas(gcf,figfile4,'epsc');
 end
 fn_fullscreen;
end

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

% NB. using fft routine as in MovingFFT

%% high pass filter by subtracting the mean;

function [Sn,ff,an,t_mean,xtra] = citeph_get_ft(time,displ)

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
 an    = 2*abs(aneg);       % nb. aneg = conj(apos) for a real signal
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