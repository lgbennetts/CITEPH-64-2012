% function Main_Attn
%
% DESCRIPTION: Compare model & data attenuations
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
% DO_FPLT  = plots on/off for each frequency
% DEL      = delete temporary files after use on/off
% 
% Model:
%
% DO_MODEL   = analyse model on/off
% Vert_Modes = number of vertical modes used for EMM method (if used)
% model_pers = abscissa (periods) used for model
% what_model = what model to use (see Main_AttnModel)           
%
% L Bennetts Sept 2013 / Adelaide

function Main_Attn

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL

if ~exist('conc','var');     conc=39; end

if ~exist('DO_SVFG','var');  DO_SVFG=0; end
if ~exist('DO_PLOT','var');  DO_PLOT=1; end
if ~exist('DO_DISP','var');  DO_DISP=0; end
if ~exist('COMM','var');     COMM   =1; end

if DO_PLOT
 if ~exist('fig','var'); fig=1; end
 if ~exist('col','var'); 
  col=' ''r.'' , ''markersize'' , 12'; 
  col_nl=' ''ro'' , ''markersize'' , 12';
  col_model=' ''c'' ';
 end
end

if ~exist('DO_FDSP','var');  DO_FDSP=0; end

%% DATA

if ~exist('DO_DATA','var'); DO_DATA=0; end

if ~exist('file_pre','var'); file_pre = 'Temp_data/a00'; end

if ~exist('T_pers','var'); T_pers=10; end

if ~exist('data_out','var') 
 data_out.name='amp-harmo-steady-1st';
 data_out.tint='t0=t0+4*Tp; t1=t0+10*Tp;'; 
end 

if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_FPLT','var');  DO_FPLT=0; end 
%if ~exist('DO_FPLT','var');  DO_FPLT='Aspec-signal'; end

HT=fn_WhatTestData(conc,'Regular',0); ht_inds=1:length(HT); %

probes=[1,2,3,8,9]; %1:10;

%% MODEL 

if ~exist('DO_MODEL','var'); DO_MODEL=1; end

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if DO_DISP; model_pers=unique(HT(2,ht_inds)); end
if ~exist('model_pers','var'); model_pers=0.65:0.05:1.85; end

what_mod = '2d BIE'; %'Boltzmann steady'; %'2d EMM'; %'2d no long';
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% NUMERICAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_MODEL
 out=Main_AttnModels(model_pers,conc,what_mod,Vert_Modes,0,0);
 count_mods=1;
 for loop_mods=1:length(out)
  if strcmp(out(loop_mods).name,'2d no long')
   trans_model(count_mods,:)=out(loop_mods).value;
   count_mods=count_mods+1;
  elseif strcmp(out(loop_mods).name,'2d long')
   trans_model(count_mods,:)=out(loop_mods).value;
   count_mods=count_mods+1;
  elseif strcmp(out(loop_mods).name,'2d BIE')
   trans_model_bie=out(loop_mods).value;
   %count_mods=count_mods+1;
  elseif strcmp(out(loop_mods).name,'periods')
   model_pers_bie(count_mods,:)=out(loop_mods).value;
  elseif strcmp(out(loop_mods).name,'Boltzmann steady')
   trans_model(count_mods,:)=out(loop_mods).value;
   count_mods=count_mods+1; 
  end
 end
 clear count_mods
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_DATA

%% INC WAVE

if COMM; cprintf(0.4*[1,1,1],'>>> incident waves\n'); end

for loop_ht=ht_inds
 if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
   '; Hs=' num2str(HT(1,loop_ht)) '\n']); end
 for probe=probes
  Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,T_pers,data_out,...
   loop_ht,probe,'LHS',DO_FPLT,'none',DO_FDSP,file_pre);
  if DO_FPLT
   pause
   close all
  end
 end
end

%% TRANSMITTED WAVE

if COMM; cprintf(0.4*[1,1,1],'>>> transmitted waves\n'); end

for loop_ht=ht_inds
 if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
   '; Hs=' num2str(HT(1,loop_ht)) '\n']); end
 for probe=probes
  Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,T_pers,data_out,...
   loop_ht,probe,'RHS',DO_FPLT,'none',DO_FDSP,file_pre);
  if DO_FPLT
   pause
   close all
  end
 end
end

%% TRANSMISSION COEFFICIENT

%%% Calculate incident amplitude as average over individual probes

inc_amp=zeros(length(dir([file_pre '*'])),1);
count_ht=1;
for loop_ht=ht_inds
 file_nms=dir([file_pre sprintf('%03g',loop_ht) '*']);
 for loop_nms=1:length(file_nms)
  file_nm=which(file_nms(loop_nms).name);
  H_vec(count_ht)=HT(1,loop_ht);
  T_vec(count_ht)=HT(2,loop_ht);
  eval(['load ' file_nm ' outputs;'])
  count_pb=1;
  for probe=probes
   out=outputs{count_pb};
   inc_amp(count_ht)=inc_amp(count_ht)+out.value;
   clear out
   count_pb=count_pb+1;
  end
  clear outputs count_pb
  count_ht=count_ht+1;
 end
end
inc_amp = inc_amp/length(probes);

%%% Calculate transmitted amplitude as average over individual probes

trans_amp=0*inc_amp;
count_ht=1;
for loop_ht=ht_inds
 file_nms=dir([file_pre sprintf('%03g',loop_ht) '*']);
 for loop_nms=1:length(file_nms)
  file_nm=which(file_nms(loop_nms).name);
  eval(['load ' file_nm ' outputs;'])
  count_pb=1;
  for probe=probes
   out=outputs{count_pb+length(probes)};
   trans_amp(count_ht)=trans_amp(count_ht)+out.value;
   clear out
   count_pb=count_pb+1;
  end
  clear outputs count_pb
  count_ht=count_ht+1;
 end
end
trans_amp = trans_amp/length(probes);

%%% Calculate ratio trans/inc

trans_coeff = trans_amp./inc_amp;

%% DELETE

if DEL
 for loop_ht=ht_inds
  file_nms=dir([file_pre sprintf('%03g',loop_ht) '*']);
  for loop_nms=1:length(file_nms)
   file_nm=which(file_nms(loop_nms).name);
   eval(['delete ' file_nm])
  end
 end
end

end % IF DO_DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%  PLOT & DISPLAY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT

if DO_PLOT
 if fig==0; fig = fn_getfig; end
 figure(fig)
 if DO_DATA
  %%% Higest amplitude waves
  [~,IA] = unique(T_vec); IA = [0,IA];
  jj_nl=[]; ct=1;
  for loop=find(diff(IA)>1)
   dum_inds=IA(loop)+1:IA(loop+1);
   if unique(H_vec(dum_inds))>1
    [~,jj_nl(ct)]=max(H_vec(dum_inds));
    jj_nl(ct)=dum_inds(jj_nl(ct)); ct=ct+1;
   end 
  end
  clear ct IA dum_inds
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  eval(['plot(T_vec,trans_coeff,' col ')'])
  if ~isempty(jj_nl)
   hold on
   eval(['plot(T_vec(jj_nl),trans_coeff(jj_nl),' col_nl ')']) 
  end % END ~ISEMPTY(jj_nl)
  prb_str=int2str(probes(1));
  for loop=2:length(probes)
   prb_str=[prb_str ', ' int2str(probes(loop))];
  end
  title(['Moving FFT:' int2str(T_pers) ' periods; interval ' ...
   data_out.tint ' probes ' prb_str],'fontsize',14)
 end
 xlabel('period [s]','fontsize',14)
 ylabel('Trans coeff','fontsize',14)
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
  saveas(figure(fig),['Figs/Attn' num2str(conc) '_' dt '@' tm '.fig']);
 end
 fn_fullscreen(fig)
end % END DO_PLOT

%% DISPLAY

if DO_DISP
 if DO_DATA
  cprintf(0.4*[1,1,1],['Periods & amplitudes:\n'])
  disp(HT(2,ht_inds))
  disp(HT(1,ht_inds)/2)
  cprintf(0.4*[1,1,1],['Transmission coefficients:\n'])
  disp(trans_coeff)
 end
 if DO_MODEL
  cprintf(0.4*[1,1,1],['Model:\n'])
  disp(model_pers)
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