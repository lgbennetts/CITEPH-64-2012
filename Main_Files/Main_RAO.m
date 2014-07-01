% function Main_RAO
%
% DESCRIPTION: Compare model & data RAOs for a single disk
%
% INPUTS:
%
% General:
%
% DO_FDSP   = comments on/off for each frequency
% DO_PLOT   = plot data on/off
% DO_STR    = strings (title/labels) for figures on/off
% DO_DISP   = print data on/off
% DO_SVFG   = save figure on/off
% PT_TYP    = 0 (angle/inc ht); 1 (gradient/k0)  
% xx        = abscissa for plot (=2 period;=3 wavelength)
% fig       = figure handle (fig=0 picks the next available figure handle)
% col       = colout/linestyle for data
% col_model = as above for model data
%
% Data:
%
% rbms     = rigid body modes: heave, surge, sway, pitch & roll
% DO_DATA  = analyse experimental data on/off
% file_pre = prefix for data files (temporary only)
% T_pers   = number of periods used for moving Fourier window
% data_out = what data are we after? 
%            data_out.name usually 'amp-harmo-steady-1st'
%            data_out.tint = steady state time interval
%            see MovingFFT.m for full description
% DO_FPLT  = plots on/off for each frequency
% DEL      = delete temporary files after use on/off
% errbars  = plot error bars (1=on, 0=off)

% 
% Model:
%
% DO_MODEL   = analyse model on/off
% Vert_Modes = number of vertical modes used for EMM method (if used)
% model_pers = abscissa (periods) used for model
% WHAT_MODEL = what model to use 
%              2d     (boundary integral eqn method)
%              2d-EMM (eigenfunction matching method)
%              3d     (eigenfunction matching method)
%
% L Bennetts Aug 2013 / Adelaide

function Main_RAO

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

if ~or(strcmp(getenv('LOGNAME'),'a1612881'),...
  strcmp(getenv('LOGNAME'),'lbennetts'))
 OTHER_USR=1;
end

%% GENERAL

if ~exist('DO_PLOT','var');  DO_PLOT=1; end
if ~exist('DO_STR','var');   DO_STR=0; end
if ~exist('DO_DISP','var');  DO_DISP=0; end
if ~exist('DO_SVFG','var');  DO_SVFG=0; end

if ~exist('PT_TYP','var');   PT_TYP =0; end

if DO_PLOT
 if ~exist('col','var'); 
  col=' ''k.'' , ''markersize'' , 12';
  col_nl=' ''ko'' , ''markersize'' , 12';
  col_model={' ''k.-'' ',' ''b'' '};
 end
end

if ~exist('xx','var'); xx=2; end

if DO_PLOT; fig = 1:3; end

if ~exist('DO_FDSP','var');  DO_FDSP=0; end

%% DATA

if ~exist('DO_DATA','var'); DO_DATA=0; end

if ~exist('file_pre','var'); file_pre = 'Temp_data/b00'; end

%if ~exist('T_pers','var'); T_pers=15; end

if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

% if ~exist('data_out','var') 
%  data_out.name='amp-harmo-steady-1st';
%  data_out.tint='t0=t0+4*Tp; t1=t0+10*Tp;'; 
% end 

if ~exist('data_out','var')
 data_out.name={'amp-harmo-steady-1st'};
 data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
end

if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_FPLT','var');  DO_FPLT=0; end 
%if ~exist('DO_FPLT','var');  DO_FPLT='Aspec-signal'; end
%if ~exist('DO_FPLT','var');  DO_FPLT='Aspec'; end

HT=fn_WhatTestData(1,'Regular',0); ht_inds=1:4; %length(HT); %length(HT); %

probes=1:10; %[1,2,3,8,9]; %10; %
rbms={'heave'}; %,'roll+pitch','surge+sway'}; %
%rbms={'heave','pitch','surge'};
%rbms={'pitch'};

if ~exist('errbars','var');  errbars=1; end

%% MODEL 

if ~exist('DO_MODEL','var');   DO_MODEL  =1; end
if ~exist('WHAT_MODEL','var'); WHAT_MODEL='2d'; end

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if DO_DISP; model_pers=HT(2,ht_inds); end
%if ~exist('model_pers','var'); model_pers=0.6:0.05:2; end
if ~exist('model_pers','var'); model_pers=fn_WhatTestData(79,'Regular',0); 
 model_pers = unique(model_pers(2,:)); end

count_rbm=1;
for loop_rbm=1:length(rbms)
 if strfind(rbms{loop_rbm},'heave')
  model_outs{count_rbm} = 'RAO-heave';
 elseif strfind(rbms{loop_rbm},'surge')
  model_outs{count_rbm} = 'RAO-surge';
 elseif strfind(rbms{loop_rbm},'pitch')
  if ~PT_TYP; model_outs{count_rbm} = 'ANG-pitch'; else
   model_outs{count_rbm} = 'RAO-pitch/k'; end
 end
 count_rbm=count_rbm+1;
end
clear count_rbm
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% NUMERICAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_MODEL
 
 %% 2d model:
 %%% Using EMM during runtime:
 
 if strfind(WHAT_MODEL,'2d-EMM')
  %cprintf('green','More checks of EMM required ?!?\n')
  if ~exist('rigid','var'); rigid=10; end
  if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid);
   Param = ModParam_def(Param,1,Vert_Modes,0,0); end
  for loop_p=1:length(model_pers)
   out = fn_2dFloe('freq',1/model_pers(loop_p),Param,'heave pitch',0,0);
   for loop_out=1:length(model_outs)
    count=1;
    while and(count>0,count<=length(out))
     if strfind(model_outs{loop_out},out(count).name)
      RAO_model_2d(loop_out,loop_p)=out(count).value;
      if strcmp(model_outs{loop_out},'ANG-pitch')
       RAO_model_2d(loop_out,loop_p)=180*RAO_model_2d(loop_out,loop_p)/pi;
      end
      count=0;
     else
      count=count+1;
     end % END IF
    end % END WHILE
    if count>length(out)
     RAO_model_2d(loop_out,loop_p)=nan;
    end
   end % END LOOP_OUT
   clear out count
  end % END LOOP P
  
  model_pers_2d=model_pers;
  
 %%% Using BIE Method & saved to file:
  
 elseif strfind(WHAT_MODEL,'2d')
  
  if exist('OTHER_USR','var')
   cprintf('red','>>> Option WHAT_MODEL=2d is not available!!! \n')
   return
  end
  
  if strcmp(getenv('LOGNAME'),'a1612881')
   basedir='/Volumes/scratch/Data/CITEPH-64-2012/Model/SingleFloe/';
  end
  disp('2d model:')
  dum_dirs = fn_FolderNames(basedir);
  for loop_d=1:length(dum_dirs)
   disp([int2str(loop_d) '. ' dum_dirs{loop_d}])
  end
  ind = str2num(input('Choose folder number... ', 's'));
  basedir=[basedir dum_dirs{ind} '/'];
  dum_dirs = fn_FileNames(basedir);
  for loop_d=1:length(dum_dirs)
   disp([int2str(loop_d) '. ' dum_dirs{loop_d}])
  end
  cprintf(0.4*[1,1,1],'xxxxxxxxxxx-Xint-Xmodes-Nfac\n')
  ind = str2num(input('Choose file number... ', 's'));
  basedir=[basedir dum_dirs{ind}];
  load(basedir,'out','for_typ','for_vec','v_info')
  cprintf(0.4*[1,1,1],[v_info{1}(1:end-1) '\n'])
  clear ind basedir dum_dirs
  if strcmp(for_typ,'period')
   model_pers_2d=for_vec; clear for_vec for_typ
  else
   cprintf('green','check the forcing type\n')
  end
  count=1;
  while and(count>0,count<=length(out))
   if strfind('wavenumber',out(count).name)
    k0=out(count).value;
    count=0;
   else
    count=count+1;
   end
  end
  for loop_out=1:length(model_outs)
   count=1;
   while and(count>0,count<=length(out))
    if strfind(model_outs{loop_out},out(count).name)
     RAO_model_2d(loop_out,:)=out(count).value;
     if strcmp(model_outs{loop_out},'ANG-pitch')
      RAO_model_2d(loop_out,:)=180*RAO_model_2d(loop_out,:)/pi;
      % Temp fix...
      if strfind(v_info{1},'1.0')
       RAO_model_2d(loop_out,:)=RAO_model_2d(loop_out,:)/.495;
      end
     elseif strcmp(model_outs{loop_out},'RAO-pitch/k')
      RAO_model_2d(loop_out,:)=RAO_model_2d(loop_out,:)./k0;
     end
     count=0;
    else
     count=count+1;
    end % END IF
   end % END WHILE
   if count>length(out)
    RAO_model_2d(loop_out,:)=nan;
   end
  end % END LOOP_OUT
  clear out count
  
 end % END IF WHAT MODEL
 
 %% 3d disk model:
 if strfind(WHAT_MODEL,'3d')
  RAO_model = Main_SingleDiskModel('freq',1./model_pers,...
   Vert_Modes,model_outs,0,DO_FDSP);
 end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_DATA
 
 if ~exist('OTHER_USR','var')
  
  if DEL
   for loop_ht=ht_inds
    file_nms=dir(['Temp_data/*']); file_nms={file_nms.name};
    file_nms(ismember(file_nms,{'.','..'})) = [];
    for loop_nms=1:length(file_nms)
     file_nm=which(file_nms{loop_nms});
     eval(['delete ' file_nm])
    end
   end
  end
  
  %% INC WAVE
  
  for loop_ht=ht_inds
   for probe=probes
    Main_SingleDiskData(HT(2,loop_ht),HT(1,loop_ht),Tpers,data_out,...
     loop_ht,probe,'WP',DO_FPLT,1,DO_FDSP,file_pre)
    if DO_FPLT
     pause
     close all
    end
   end
  end
  
  %% RIGID BODY MODES
  
  for loop_ht=ht_inds
   for loop_rbm=1:length(rbms)
    if strfind(rbms{loop_rbm},'+')
     ind = strfind(rbms{loop_rbm},'+');
     probe=rbms{loop_rbm}(1:ind-1);
     Main_SingleDiskData(HT(2,loop_ht),HT(1,loop_ht),T_pers,data_out,...
      loop_ht,probe,'RBM',DO_FPLT,'none',DO_FDSP,file_pre)
     probe=rbms{loop_rbm}(ind+1:end);
     Main_SingleDiskData(HT(2,loop_ht),HT(1,loop_ht),T_pers,data_out,...
      loop_ht,probe,'RBM',DO_FPLT,'none',DO_FDSP,file_pre)
    else
     probe=rbms{loop_rbm};
     Main_SingleDiskData(HT(2,loop_ht),HT(1,loop_ht),Tpers,data_out,...
      loop_ht,probe,'RBM',DO_FPLT,'none',DO_FDSP,file_pre)
    end % END IF '+'
    if DO_FPLT
     pause
     close all
    end
   end
  end
  
  %% RESPONSE AMPLITUDE OPERATORS (RAOs)
  
  %%% Calculate incident amplitude as average over individual probes
  
  if isempty(probes)
   inc_amp=ones(length(ht_inds),1);
  else
   inc_amp=zeros(length(ht_inds),1);
   count_ht=1;
   for loop_ht=ht_inds
    file_nm=fn_get_filenm(loop_ht,file_pre);
    eval(['load ' file_nm ' outputs;'])
    count_pb=1;
    for probe=probes
     out=outputs{count_pb};
     amps_hld(count_pb)=out.value;
     %inc_amp(count_ht)=inc_amp(count_ht)+out.value;
     clear out
     count_pb=count_pb+1;
    end % END PROBES
    inc_amp(count_ht)=mean(amps_hld);
    inc_std(count_ht)=std(amps_hld);
    clear outputs count_pb amps_hld
    count_ht=count_ht+1;
   end % END LOOP HT_INDS
   %inc_amp = inc_amp/length(probes);
  end % END IF ISEMPTY(PROBES)
  
  %%% RAOs
  
  count_ht=1;
  for loop_ht=ht_inds
   file_nm=fn_get_filenm(loop_ht,file_pre);
   eval(['load ' file_nm ' outputs;'])
   count_rbm=1;
   for loop_rbm=1:length(rbms)
    probe=rbms{loop_rbm};
    if strfind(probe,'+')
     out=outputs{length(probes)+count_rbm};
     dum_RAO=out.value;
     count_rbm=count_rbm+1;
     out=outputs{length(probes)+count_rbm};
     RAO(count_ht,loop_rbm)=abs(out.value + 1i*dum_RAO);
     count_rbm=count_rbm+1; clear dum_RAO
    else
     probe=rbms{loop_rbm};
     out=outputs{length(probes)+count_rbm};
     RAO(count_ht,loop_rbm)=out.value;
     count_rbm=count_rbm+1;
    end % END IF '+'
    if PT_TYP
     if or(or(~isempty(strfind(probe,'surge')),...
       ~isempty(strfind(probe,'sway'))),...
       ~isempty(strfind(probe,'heave')))
      dumRAO_std(1,count_ht,loop_rbm) = (RAO(count_ht,loop_rbm)+0.001)./...
       (inc_amp(count_ht)-0.001);
      dumRAO_std(2,count_ht,loop_rbm) = (RAO(count_ht,loop_rbm)-0.001)./...
       (inc_amp(count_ht)+0.001);
      RAO(count_ht,loop_rbm)       = RAO(count_ht,loop_rbm)/...
       inc_amp(count_ht);
     else
      k0=2*pi/HT(3,ht_inds(count_ht))/.99;
      dumRAO_std(1,count_ht,loop_rbm) = tan(pi*(RAO(count_ht,loop_rbm)+0.1)/180)./...
       (inc_amp(count_ht)-0.001)/k0;
      dumRAO_std(2,count_ht,loop_rbm) = tan(pi*(RAO(count_ht,loop_rbm)-0.1)/180)./...
       (inc_amp(count_ht)+0.001)/k0;
      RAO(count_ht,loop_rbm)=tan(pi*RAO(count_ht,loop_rbm)/180)/...
       inc_amp(count_ht)/k0; clear k0
     end % END IF HEAVE OR SURGE OR SWAY
    else
     dumRAO_std(1,count_ht,loop_rbm) = RAO(count_ht,loop_rbm)./(inc_amp-inc_std);
     dumRAO_std(2,count_ht,loop_rbm) = RAO(count_ht,loop_rbm)./(inc_amp+inc_std);
     RAO(count_ht,loop_rbm)=RAO(count_ht,loop_rbm)/inc_amp(count_ht);
    end % END IF PT_TYP
    clear out
   end
   clear outputs
   count_ht=count_ht+1;
  end
  
  %% DELETE
  
  if DEL
   for loop_ht=ht_inds
    file_nm=fn_get_filenm(loop_ht,file_pre);
    eval(['delete ' file_nm '.mat'])
    %%
   end
  end
  
 %%%%%%%%%%%%%%%%%%%%%
 %%%  IF OTHER_USR 
 %%%%%%%%%%%%%%%%%%%%%
 else 
  
  cprintf('red','>>> DO_DATA: recovering data points only \n')
  cprintf('red','    i.e. no calculations \n')
  
  dum_RAO = zeros(7,length(rbms));
  dum_RAO_std = zeros(2,7,length(rbms));
  
  for lp=1:length(rbms)
   if ~isempty(strfind(rbms{lp},'heave'))
    dum_dat=load('Data/RAO_Heave');
    dum_RAO(:,loop_rbm)=dum_dat.Hv;
    dum_RAO_std(:,:,loop_rbm)=dum_dat.Hv_std;
   elseif ~isempty(strfind(rbms{lp},'pitch'))
    dum_dat=load('Data/RAO_Pitch');
    dum_RAO(:,loop_rbm)=dum_dat.Pt;
    dum_RAO_std(:,:,loop_rbm)=dum_dat.Pt_std;
   elseif ~isempty(strfind(rbms{lp},'surge'))
    dum_dat=load('Data/RAO_Surge');
    dum_RAO(:,loop_rbm)=dum_dat.Sg;
    dum_RAO_std(:,:,loop_rbm)=dum_dat.Sg_std;
   end
  end
  T_all = dum_dat.T; H_all = dum_dat.H;
  ct=1;
  for lp=1:length(ht_inds)
   inds=find(and(T_all==HT(2,ht_inds(lp)),H_all==HT(1,ht_inds(lp))));
   for lpp=1:length(inds)
    RAO(ct,:)=dum_RAO(inds(lpp),:);
    RAO_std(:,ct,:)=dum_RAO_std(:,inds(lpp),:);
    ct=ct+1;
   end
  end
  clear dum_dat H_all T_all dum_RAO dum_RAO_std
  
 end % IF EXIST OTHER_USR
 
end % IF DO_DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%  PLOT & DISPLAY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT

if DO_PLOT
 if ~exist('fig','var'); fig = fn_getfig(length(rbms)); end
 %%% Higest amplitude waves
 [~,IA] = unique(HT(2,ht_inds)); IA = [0,IA];
 jj_nl=[]; ct=1;
 for loop=find(diff(IA)>1)
  dum_inds=IA(loop)+1:IA(loop+1);
  if unique(HT(1,dum_inds))>1
   [~,jj_nl(ct)]=max(HT(1,dum_inds));
   jj_nl(ct)=dum_inds(jj_nl(ct)); ct=ct+1;
  end
 end
 clear ct IA dum_inds
 %%%%%%%%%%%%%%%%%%%%%%%%%%%
 for loop_rbm=1:length(rbms)
  figure(fig(loop_rbm)); hold on
  %%%%%%%%%%%%
  %%% DATA %%%
  %%%%%%%%%%%%
  if DO_DATA
   eval(['plot(HT(xx,ht_inds),RAO(:,loop_rbm),' col ')'])
   if ~isempty(jj_nl)
    %hold on
    eval(['plot(HT(xx,ht_inds(jj_nl)),RAO(jj_nl,loop_rbm),' col_nl ')'])
    %hold off
   end
   if errbars
    fn_bar(HT(xx,ht_inds(setdiff(1:length(ht_inds),jj_nl))),...
     squeeze(RAO_std(:,setdiff(1:length(ht_inds),jj_nl),loop_rbm)),0.01,1);
    fn_bar(HT(xx,ht_inds(jj_nl)),...
     squeeze(RAO_std(:,jj_nl,loop_rbm)),0.01,0);
   end
   if DO_STR
    prb_str=int2str(probes(1));
    for loop=2:length(probes)
     prb_str=[prb_str ', ' int2str(probes(loop))];
    end
    title([rbms{loop_rbm} ': ' Tpers ' periods; interval ' ...
     data_out.tint ' probes ' prb_str],'fontsize',14)
   end % END DO_STR
  end
  if DO_STR
   if xx==2
    xlabel('period [s]','fontsize',14)
   elseif xx==3
    xlabel('wavelength/floe diameter','fontsize',14)
   end
   ylabel(['RAO-' rbms{loop_rbm}],'fontsize',14)
  end %END DO_STR
  %%%%%%%%%%%%%
  %%% MODEL %%%
  %%%%%%%%%%%%%
  if DO_MODEL
   hold on
   if strfind(WHAT_MODEL,'2d')
    if xx==2
     eval(['plot(model_pers_2d,RAO_model_2d(loop_rbm,:),' col_model{1} ')'])
    elseif xx==3
     eval(['plot(2*pi./k0/.99,RAO_model_2d(loop_rbm,:),' col_model{1} ')'])
    end
   end
   if strfind(WHAT_MODEL,'3d')
    eval(['plot(model_pers,RAO_model(loop_rbm,:),' col_model{2} ')'])
   end
   hold off
  end % END IF DO_MODEL
  if DO_SVFG
   clockout=clock;
   dt=[sprintf('%0.4d',clockout(1)) '-' sprintf('%0.2d',clockout(2)) '-' sprintf('%0.2d',clockout(3))];
   tm=[sprintf('%0.2d',clockout(4)) ':' sprintf('%0.2d',clockout(5)) ':' sprintf('%0.2d',floor(clockout(6)))];
   saveas(figure(fig(loop_rbm)),['Figs/RAO_' rbms{loop_rbm} '_' dt '@' tm '.fig']);
  end
  fn_fullscreen(fig(loop_rbm))
 end % END LOOP_RBM
end % END DO_PLOT

%% DISPLAY

if DO_DISP
 cprintf(0.4*[1,1,1],['Periods:\n'])
 disp(HT(2,ht_inds))
 for loop_rbm=1:length(rbms)
  if and(~exist('k0','var'),...
    or(or(or(~isempty(strfind(rbms{loop_rbm},'surge')),...
     ~isempty(strfind(rbms{loop_rbm},'sway'))),...
     ~isempty(strfind(rbms{loop_rbm},'pitch'))),...
     ~isempty(strfind(rbms{loop_rbm},'roll'))))
   bed=3.1; k0=zeros(1,length(ht_inds));
   count_ht=1;
   for loop_ht=ht_inds
    sig = ((2*pi/HT(2,loop_ht))^2)/9.81;
    k0(count_ht) = CalcRealRoot_PWC([bed, sig], ...
     'FS_DispRel_PWC', 'UppLimReal_FS_PWC', 1e-16);
    count_ht=count_ht+1;
   end
   clear count_ht sig
  end % end calc k0
 end
 
 for loop_rbm=1:length(rbms)
  count_ht=1;
  for loop_ht=ht_inds
   if or(~isempty(strfind(rbms{loop_rbm},'surge')),...
     ~isempty(strfind(rbms{loop_rbm},'sway')))
   disp_vec(:,count_ht)=[RAO(count_ht,loop_rbm);...
    coth(k0(count_ht)*bed)];
  elseif or(~isempty(strfind(rbms{loop_rbm},'pitch')),...
    ~isempty(strfind(rbms{loop_rbm},'roll')))
   disp_vec(:,count_ht)=[RAO(count_ht,loop_rbm)*inc_amp(count_ht);...
    tan(pi*RAO(count_ht,loop_rbm)*inc_amp(count_ht)/180)/k0(count_ht);...
    RAO(count_ht,loop_rbm)];
  else
   disp_vec(:,count_ht)=[RAO(count_ht,loop_rbm)];
  end
  count_ht=count_ht+1;
  end % end loop_ht
  if or(~isempty(strfind(rbms{loop_rbm},'surge')),...
     ~isempty(strfind(rbms{loop_rbm},'sway')))
   cprintf(0.4*[1,1,1],[rbms{loop_rbm} ' RAO, eccentricity:\n'])
  elseif or(~isempty(strfind(rbms{loop_rbm},'pitch')),...
    ~isempty(strfind(rbms{loop_rbm},'roll')))
   cprintf(0.4*[1,1,1],[rbms{loop_rbm} ' angle (degs), tan(ang)/k, angle/amp:\n'])
  else
   cprintf(0.4*[1,1,1],[rbms{loop_rbm} ' RAO:\n'])
  end
  disp(disp_vec)
  if DO_MODEL
   cprintf(0.4*[1,1,1],['Model ' model_outs{loop_rbm} ':\n'])
   disp(RAO_model(loop_rbm,:))
  end % END IF DO_MODEL
  clear disp_vec
 end % for loop_rbm
end % end DO_DISP

return

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function file_nm=fn_get_filenm(run,file_pre)

if ~exist('file_pre','var'); file_pre = 'Temp_data/s13'; end

if run > 99
  file_nm = [file_pre int2str(run)];
 elseif run > 9
  file_nm = [file_pre '0' int2str(run)];
 else
  file_nm = [file_pre '00' int2str(run)];
end

return 