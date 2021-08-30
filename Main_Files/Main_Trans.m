% function Main_Trans
%
% DESCRIPTION: Compare model & data wave transmissions
%
% INPUTS:
%
% General:
%
% conc      = either 39 or 79
% DO_FDSP   = comments on/off for each frequency
% DO_PLOT   = plot data on/off
% DO_DISP   = print data on/off
% DO_SVFG   = save figure on/off
% COMM      = comments on/off
% EGY       = transmitted energy rather than amplitude (1=on, 0=off)
% fig       = figure handle (fig=0 picks the next available figure handle)
% col       = colout/linestyle for data
% col_model = as above for model data
%
% Data:
%
% probes   = probes used for analysis
% DO_DATA  = analyse experimental data on/off
% file_pre = prefix for data files (temporary only)
% T_pers   = number of periods used for moving Fourier window
% data_out = what data are we after?
%            data_out.name usually 'amp-harmo-steady-1st'
%            data_out.tint = steady state time interval
%            see MovingFFT.m for full description
% t_meth   = data to compare transmitted waves to
%            'inc' using incident wave
%            'calib' using calibration test
% errbars  = plot error bars (1=on, 0=off)
% DO_FPLT  = plots on/off for each frequency
% DEL      = delete temporary files after use on/off
%
% Model:
%
% DO_MODEL   = analyse model on/off
% Vert_Modes = number of vertical modes used for EMM method (if used)
% model_pers = abscissa (periods) used for model
% what_mod   = what model to use (see Main_AttnModel)
%
% L Bennetts Sept 2013 / Adelaide

function Main_Trans(my_inputs)

if exist('my_inputs','var'); eval(my_inputs); end

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

if ~(strcmp(getenv('LOGNAME'),'a1612881')||...
  strcmp(getenv('LOGNAME'),'lbennetts')||...
  strcmp(getenv('LOGNAME'),'a1229158'))
 OTHER_USR=1;
end

%% GENERAL

if ~exist('conc','var') || isempty(conc);     conc=39; end
if ~exist('WaveType','var');     WaveType= 'Regular';end %'Irregular %'Regular'

if ~exist('DO_SVFG','var');  DO_SVFG=0; end
if ~exist('DO_PLOT','var');  DO_PLOT=1; end
if ~exist('DO_DISP','var');  DO_DISP=0; end
if ~exist('COMM','var');     COMM   =1; end
if ~exist('EGY','var');      EGY    =1; end
if ~exist('DO_STR','var');   DO_STR =1; end


if ~exist('PLOT_TYP','var'); PLOT_TYP ='trans_coeff'; end

if DO_PLOT
 if ~exist('fig','var'); fig=1; end
 if ~exist('col','var');
%   col=' ''k.'' , ''markersize'' , 12';
%   col_nl=' ''ko'' , ''markersize'' , 12';
%   col_ii=' ''kx'' , ''markersize'' , 12';
%   col_nl_ii=' ''ko'' , ''markersize'' , 12';
%   coll=' ''b.'' , ''markersize'' , 12';
%   collnl=' ''bo'' , ''markersize'' , 12';
%   col_model=' ''b-.'' ';
  col=' ''k.'' , ''markersize'' , 12';
  col_nl=' ''ko'' , ''markersize'' , 12';
  col_ii=' ''kx'' , ''markersize'' , 12';
  col_nl_ii=' ''ko'' , ''markersize'' , 12';
  coll=' ''k.'' , ''markersize'' , 12';
  collnl=' ''ko'' , ''markersize'' , 12';
  col_model=' ''k-.'' ';
%   col=' ''r.'' , ''markersize'' , 12';
%   col_nl=' ''ro'' , ''markersize'' , 12';
%   col_ii=' ''rx'' , ''markersize'' , 12';
%   col_nl_ii=' ''ro'' , ''markersize'' , 12';
%   coll=' ''r.'' , ''markersize'' , 12';
%   collnl=' ''ro'' , ''markersize'' , 12';
%   col_model=' ''r-.'' ';
%col_model=' ''go'' , ''markersize'' , 12';
 end
end

%% DATA

if ~exist('DO_DATA','var'); DO_DATA=1; end

if ~exist('t_meth','var');  t_meth='calib'; end %'inc'; 'calib'

if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end %'Tpers=10;'; end %

if ~exist('data_out','var')
 data_out.name={'amp-harmo-steady-1st'};
 %data_out.name={'amp-harmo-steady-1st','amp-harmo-steady-2nd'};
 %data_out.tint='t0=t0+4*Tp; t1=t0+10*Tp;'; %t1=min(t0+10*Tp,min(t_vec));'; % 
 %data_out.tint='t0=t0+(T_pers/2)*Tp; t1=min(t_vec)-(T_pers/2)*Tp;'; 
 data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
end

if ~exist('DEL','var');      DEL=0; end

if ~exist('DO_FDSP','var');  DO_FDSP=0; end
if ~exist('DO_FPLT','var');  DO_FPLT=0; end
%if ~exist('DO_FPLT','var');  DO_FPLT='Aspec-signal'; end 
%if ~exist('DO_FPLT','var');  DO_FPLT='vid'; end %

if ~exist('errbars','var');  errbars=1; end

HT =fn_WhatTestData(conc,WaveType,0); ht_inds=1:length(HT(1,:)); %1; %
% if s trfind(t_meth,'calib')  
%  HT0=fn_WhatTestData(conc,['Regular-' t_meth],0);
%  intersect(intersect(HT(1,:),HT0(1,:)),intersect(HT(2,:),HT0(2,:)))
% end

probes=1:10; %1; % [1:3,7,10]; %[1,2,3,8,9,10]; % 

if ~exist('file_pre','var'); file_pre = 'Temp_data/a00'; end

%% MODEL

if ~exist('DO_MODEL','var'); DO_MODEL=0; end

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if DO_DISP; model_pers=0.3:0.02:2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end %model_pers=unique(HT(2,ht_inds)); end % 0.6:0.02:2; end %

what_mod = 'Boltzmann steady';%'2d EMM';%'Boltzmann steady';%'2d EMM' ;%'Boltzmann steady';%'2d EMM'; %'Boltzmann steady'; %'2d EMM'; %'Boltzmann steady' ;%'2d EMM';
RemScat=0;

% Consistency for other users:

if exist('OTHER_USR','var')
 if ~strcmp(t_meth,'calib')
  cprintf('r','>>> Setting t_meth=calib\n')
  t_meth='calib';
 end
 if EGY~=1
  cprintf('r','>>> Setting EGY=1\n')
  EGY    =1;
 end
 if ~strcmp(Tpers,'Tpers=fn_Tpers(Tp);')
  cprintf('r','>>> Setting Tpers=fn_Tpers(Tp)\n')
  Tpers='Tpers=fn_Tpers(Tp);';
 end
 if ~strcmp(data_out.tint,'[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);')
  cprintf('r','>>> Setting data_out.tint=[t0,t1] = fn_tint(tvec,t0,Tp,Tpers)\n')
  data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
 end
 if length(probes)~=10
  cprintf('r','>>> Setting probes=1:10\n')
  probes=1:10;
 end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if length(data_out.name)==1

    cprintf('red','>>> DO_DATA: recovering data points only \n')
    cprintf('red','    i.e. no calculations \n')

    if ~strcmp(PLOT_TYP,'trans_coeff')
    cprintf('red','>>> transmission coefficients only!!! \n')
    cprintf('red','    set PLOT_TYP=''trans_coeff'' \n')
    return
    end
    
    DirBase = 'Data/Gen/';
    File_NM = strcat(DirBase,'Trans',int2str(conc),WaveType(1:3));
    
    if (~exist(strcat(File_NM,'.mat'), 'file'))
        fn_GenerateTransMatrices(conc,WaveType,t_meth)
    end
    dum_dat = load(File_NM);

%     if conc==39
%     dum_dat = load('Data/Trans38');
%     elseif conc==79
%     dum_dat = load('Data/Trans77');
%     else
%     cprintf('red','>>> need conc=39 or 79 \n')
%     return
%     end 

    dum_trans     = dum_dat.trans;
    dum_trans_err = dum_dat.trans_std;
    T_all = dum_dat.T; H_all = dum_dat.H;
    
    %Remove the likely errors
%     trans_coeff = dum_trans(dum_trans < 1);
%     trans_coeff_err = dum_trans_err(:,dum_trans < 1);
%     T_vec = T_all(dum_trans < 1);
%     H_vec = H_all(dum_trans < 1);

    ct=1;
    for lp=1:length(ht_inds)
    inds=find(and(T_all==HT(2,ht_inds(lp)),H_all==HT(1,ht_inds(lp))));
    for lpp=1:length(inds)
     trans_coeff(ct,:)=dum_trans(inds(lpp));
     trans_coeff_err(:,ct,:)=dum_trans_err(:,inds(lpp));
     H_vec(ct) = H_all(inds(lpp)); T_vec(ct) = T_all(inds(lpp));
     ct=ct+1;
    end
    end

    clear dum_dat H_all T_all dum_trans dum_trans_err

    else % Higher-order harmonics

    cprintf('red','>>> Higher-order harmonics not available \n')
    return

end % END IF length(data_out.name)==1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% NUMERICAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_MODEL
 
 if and(exist('OTHER_USR','var'),strcmp(what_mod,'2d BIE'))
  cprintf('red','>>> Option WHAT_MODEL=2d BIE is not available!!! \n')
  return
 end
 
 if conc==39
  dum_c = 100*pi*(0.495^2)/2;
 elseif conc==79
  dum_c = 100*pi*(0.495^2);
 end
 out=Main_AttnModels(model_pers,dum_c,what_mod,RemScat,Vert_Modes,DO_FDSP,0);
 count_mods=1;
 for loop_mods=1:length(out)
  if strcmp(out(loop_mods).name,'2d no long')
   trans_model(count_mods,:)=out(loop_mods).value;
   if EGY; trans_model(count_mods,:)=trans_model(count_mods,:).^2; end
   count_mods=count_mods+1;
  elseif strcmp(out(loop_mods).name,'2d long')
   trans_model(count_mods,:)=out(loop_mods).value;
   if EGY; trans_model(count_mods,:)=trans_model(count_mods,:).^2; end
   count_mods=count_mods+1;
  elseif strcmp(out(loop_mods).name,'2d BIE')
   trans_model_bie=out(loop_mods).value;
   if EGY; trans_model_bie=trans_model_bie.^2; end
  elseif strcmp(out(loop_mods).name,'periods')
   model_pers_bie(count_mods,:)=out(loop_mods).value;
  elseif strcmp(out(loop_mods).name,'Boltzmann steady')
   trans_model(count_mods,:)=out(loop_mods).value;
   if EGY; trans_model(count_mods,:)=trans_model(count_mods,:).^2; end
   count_mods=count_mods+1;
  elseif strcmp(out(loop_mods).name,'Rows')
   trans_model(count_mods,:)=out(loop_mods).value;
   if EGY; trans_model(count_mods,:)=trans_model(count_mods,:).^2; end
   count_mods=count_mods+1;
  end
 end
 clear count_mods
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%  PLOT & DISPLAY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT

if DO_PLOT
 if fig==0; fig = fn_getfig; end
 figure(fig); hold on
 %%%%%%%%%%%%
 %%% DATA %%%
 %%%%%%%%%%%%
 if DO_DATA
  if ~exist('eps_vec','var'); hd(1)=subplot(1,1,1); hold on;
  else  hd(1)=subplot(2,1,1); hold on; hd(2)=subplot(2,1,2); hold on; end
  %%% Higest amplitude waves
  [~,IA] = unique(T_vec);
  if size(IA,1)==1
     IA = [0,IA];
  else
     IA = [0,IA.'];
  end
  jj_nl=[]; ct=1;
  for loop=find(diff(IA)>1)
   dum_inds=IA(loop)+1:IA(loop+1);
   if length(unique(H_vec(dum_inds)))>1
    [~,jj_nl(ct)]=max(H_vec(dum_inds));
    jj_nl(ct)=dum_inds(jj_nl(ct)); ct=ct+1;
   end
  end
  clear ct IA dum_inds
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  for lp=1:length(hd)
   if lp==1; xx=T_vec; elseif lp==2; xx=eps_vec; end
   if and(errbars,length(hd)==1)
    fn_bar(xx(setdiff(1:length(xx),jj_nl)),...
     trans_coeff_err(:,setdiff(1:length(xx),jj_nl)),0.01,1);
    fn_bar(xx(jj_nl),trans_coeff_err(:,jj_nl),0.01,0);
   end
   if strcmp(PLOT_TYP,'trans_coeff')
    eval(['p(1,lp)=plot(hd(lp),xx,trans_coeff,' col ');'])
    if exist('inc_coeff','var'); eval(['p(2,lp)=plot(hd(lp),xx,inc_coeff,' coll ');']); end
   else
    eval(['p(1,lp)=plot(hd(lp),xx,trans_std./trans_amp,' col ');'])
    eval(['p(1,lp)=plot(hd(lp),xx,inc_std./inc_amp,' col_ii ');'])
   end
   if ~isempty(jj_nl)
    if strcmp(PLOT_TYP,'trans_coeff')
     eval(['plot(hd(lp),xx(jj_nl),trans_coeff(jj_nl),' col_nl ')'])
     if exist('inc_coeff','var'); eval(['plot(hd(lp),xx(jj_nl),inc_coeff(jj_nl),' collnl ')']); end
    else
     eval(['plot(hd(lp),xx(jj_nl),trans_std(jj_nl)./trans_amp(jj_nl),' col_nl ')'])
     eval(['plot(hd(lp),xx(jj_nl),inc_std(jj_nl)./inc_amp(jj_nl),' col_nl_ii ')'])
    end
   end % END ~ISEMPTY(jj_nl)
   prb_str=int2str(probes(1));
   for loop=2:length(probes)
    prb_str=[prb_str ', ' int2str(probes(loop))];
   end
   if DO_STR
    if lp==1;     xlabel(hd(lp),'period [s]','fontsize',14);
    elseif lp==2; xlabel(hd(lp),'steepness','fontsize',14); end
   end % END DO_STR
  end
  if DO_STR
   if length(data_out.name)==1
   if strfind(data_out.name{1},'1st')
    ttl=title(['Moving FFT:' Tpers ' periods; interval ' ...
     data_out.tint ' probes ' prb_str],'fontsize',14);
    set(ttl,'interpreter','none')
   elseif strfind(data_out.name{1},'2nd')
    ttl=title(['Moving FFT: 2nd harmo, ' Tpers ' periods; interval ' ...
     data_out.tint ' probes ' prb_str],'fontsize',14);
    set(ttl,'interpreter','none')
   end
    if strcmp(PLOT_TYP,'trans_coeff')
     if ~EGY
      ylabel(hd,['Trans coeff (v ' t_meth ')'],'fontsize',14)
     else
      ylabel(hd,['Trans energy (v ' t_meth ')'],'fontsize',14)
     end
    else
     ylabel(hd,['std/amp (v ' t_meth ')'],'fontsize',14)
    end
   elseif length(data_out.name)==2
    ttl=title(hd(1),['Moving FFT:' int2str(Tpers) ' periods; interval ' ...
     data_out.tint ' probes ' prb_str],'fontsize',14);
    set(ttl,'interpreter','none')
    for lp=1:length(hd)
     ylabel(hd(lp),['A_{2}/A_{1} (v ' t_meth ')'],'fontsize',14)
    end
   end
  end % END IF DO_STR
  if length(data_out.name)==2
   legend(p(:,1),{'with floes','no floes'},'Location','NorthEast');
  end
 end
 %%%%%%%%%%%%%
 %%% MODEL %%%
 %%%%%%%%%%%%%
 if DO_MODEL
  hold on
  if strfind(what_mod,'2d BIE')
   eval(['plot(model_pers_bie,trans_model_bie,' col_model ')'])
  end
  if exist('trans_model','var')
   eval(['plot(model_pers,trans_model,' col_model ')'])
  end
  hold off
 end % END IF DO_MODEL
 if DO_SVFG
  clockout=clock;
  dt=[sprintf('%0.4d',clockout(1)) '-' sprintf('%0.2d',clockout(2)) '-' sprintf('%0.2d',clockout(3))];
  tm=[sprintf('%0.2d',clockout(4)) ':' sprintf('%0.2d',clockout(5)) ':' sprintf('%0.2d',floor(clockout(6)))];
  if length(data_out.name)==1
   saveas(figure(fig),['Figs/Trans' num2str(conc) '_' dt '@' tm '.fig']);
  elseif length(data_out.name)==2
   saveas(figure(fig),['Figs/Harmo2nd_' num2str(conc) '_' dt '@' tm '.fig']);
  end
 end
 %fn_halfscreen(fig)
end % END DO_PLOT

%% DISPLAY

if DO_DISP
 if DO_DATA
  cprintf(0.4*[1,1,1],['Data: periods & amplitudes:\n'])
  disp(T_vec) %disp(HT(2,ht_inds))
  disp(H_vec) %disp(HT(1,ht_inds)/2)
  cprintf(0.4*[1,1,1],['Transmission coefficients:\n'])
  disp(trans_coeff.')
 end
 if DO_MODEL
  cprintf(0.4*[1,1,1],['Model: periods\n'])
  disp(model_pers)
  cprintf(0.4*[1,1,1],['Transmission coefficients:\n'])
  disp(trans_model)
 end % END IF DO_MODEL
end % end DO_DISP

return

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

% function file_nm=fn_get_filenm(run,file_pre)
%
% if ~exist('file_pre','var'); file_pre = 'Temp_data/s13'; end
%
% if run > 99
%   file_nm = [file_pre int2str(run)];
%  elseif run > 9
%   file_nm = [file_pre '0' int2str(run)];
%  else
%   file_nm = [file_pre '00' int2str(run)];
% end
%
% return
