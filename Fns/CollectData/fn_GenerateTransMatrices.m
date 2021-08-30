% function fn_GenerateTransMatrices
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
%
% General:
%
% conc      = either 39 or 79
% DO_FDSP   = comments on/off for each frequency
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
%            'calib' using calibration test  (this matches paper results)
% errbars  = plot error bars (1=on, 0=off)
% DO_FPLT  = plots on/off for each frequency
% DEL      = delete temporary files after use on/off
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function fn_GenerateTransMatrices(conc,WaveType,t_meth)

if exist('my_inputs','var'); eval(my_inputs); end

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Which user?
[user_number,user_name] = citeph_get_user();

if user_number == 0
    OTHER_USR=1;
end


%% GENERAL

if ~exist('conc','var') conc=79; end
if ~exist('WaveType','var');     WaveType= 'Regular';end; %'Irregular'%'Regular'

if ~exist('DO_SVFG','var');  DO_SVFG=0; end
if ~exist('COMM','var');     COMM   =0; end
if ~exist('EGY','var');      EGY    =1; end
if ~exist('DO_STR','var');   DO_STR =1; end


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

if ~exist('DO_FDSP','var');  DO_FDSP=1; end
if ~exist('DO_FPLT','var');  DO_FPLT=0; end

HT =fn_WhatTestData(conc,WaveType,0); ht_inds=1:length(HT(1,:)); %1; %

%}{
% HT = HT(:,end);
% if s trfind(t_meth,'calib')  
%  HT0=fn_WhatTestData(conc,['Regular-' t_meth],0);
%  intersect(intersect(HT(1,:),HT0(1,:)),intersect(HT(2,:),HT0(2,:)))
% end

probes=1:10; %1; % [1:3,7,10]; %[1,2,3,8,9,10]; % 

if ~exist('file_pre','var'); file_pre = 'Temp_data/a00'; end


% Consistency for other users:

if exist('OTHER_USR','var')
  cprintf('red','>>> USER NOT DETECTED, EITHER UPDATE INFO IN citeph_get_user and citeph_user_specifics to input data \n')  
  cprintf('red','    OR RELOAD TRANSMISSION MATRICES PROVIDED ON GITHUB \n') 
  error('NO USER DEFINED / LACK OF DATA REPOSITRY INFORMATION')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%% DELETE - ALWAYS DELETE BEFORE WRITING - common error, sometimes we want the outputs from previous runs though
file_nms = dir(['Temp_data/*']);
file_nms={file_nms.name};
file_nms(ismember(file_nms,{'.','..'})) = [];
for loop_nms=1:length(file_nms)
  file_nm=strcat('Temp_data/',file_nms{loop_nms});
  eval(['delete ' file_nm])
end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% TRANSMISSION COEFFICIENT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if length(data_out.name)==1
   
   %% TRANSMITTED WAVE
   
   if strcmp(t_meth,'calib-check')
    
    if COMM; cprintf(0.4*[1,1,1],'>>> calibration waves LHS\n'); end
    
    for loop_ht=ht_inds
     if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
       '; Hs=' num2str(HT(1,loop_ht)) '\n']); end
     for probe=probes
      Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,Tpers,data_out,...
       loop_ht,probe,WaveType,'LHS-calib',DO_FPLT,'none',DO_FDSP,file_pre);
      if DO_FPLT
       pause
       close all
      end
     end
    end
    
   else
    
    if COMM; cprintf(0.4*[1,1,1],'>>> transmitted waves\n'); end
    
    for loop_ht=ht_inds
     if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
       '; Hs=' num2str(HT(1,loop_ht)) '\n']); end
     for probe=probes
      Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,Tpers,data_out,...
       loop_ht,probe,WaveType,'RHS',DO_FPLT,'none',DO_FDSP,file_pre);
      if DO_FPLT
       pause
       close all
      end
     end
    end
    
   end
   
   %% INC WAVE
   
   if strcmp(t_meth,'inc')
    
    if COMM; cprintf(0.4*[1,1,1],'>>> incident waves\n'); end
    
    for loop_ht=ht_inds
     if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
       '; Hs=' num2str(HT(1,loop_ht)) '\n']); end
     for probe=probes
      Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,Tpers,data_out,...
       loop_ht,probe,WaveType,'LHS',DO_FPLT,'none',DO_FDSP,file_pre);
      if DO_FPLT
       pause
       close all
      end
     end
    end
    
   elseif strfind(t_meth,'calib')
    
    if COMM; cprintf(0.4*[1,1,1],'>>> calibration waves\n'); end
    
    for loop_ht=ht_inds
     if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
       '; Hs=' num2str(HT(1,loop_ht)) '\n']); end
     for probe=probes
      Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,Tpers,data_out,...
       loop_ht,probe,WaveType,'RHS-calib',DO_FPLT,'none',DO_FDSP,file_pre);
      if DO_FPLT
       pause
       close all
      end
     end
    end
    
   else
    
    cprintf('green','>>> NEED VALID INPUT FOR t_meth !!!!! \n');
    
   end
   
   %% TRANSMISSION COEFFICIENT
   
   %%% Calculate transmitted amplitude as average over individual probes
   
   inc_amp=zeros(length(dir([file_pre '*'])),1);
   inc_std=inc_amp;
   trans_amp=0*inc_amp;
   trans_std=trans_amp;
   
   count_ht=1;
   for loop_ht=ht_inds
    file_nms=dir([file_pre sprintf('%03g',loop_ht) '*']);
    for loop_nms=1:length(file_nms)
     file_nm=which(file_nms(loop_nms).name);
     amps_hld       =zeros(1,length(probes));
     eval(['load ' file_nm ' outputs;']);
     %outputs = load(file_nm);
     count_pb=1;
     for probe=probes
      out=outputs{count_pb};
      amps_hld(count_pb)=out.value;
      clear out
      count_pb=count_pb+1;
     end
     trans_amp(count_ht)=mean(amps_hld);
     trans_std(count_ht)=std(amps_hld);
     clear outputs count_pb amps_hld
     count_ht=count_ht+1;
    end
   end
   
   %%% Calculate incident amplitude as average over individual probes
   
   count_ht=1;
   for loop_ht=ht_inds
    file_nms=dir([file_pre sprintf('%03g',loop_ht) '*']);
    if and(~isempty(strfind(t_meth,'calib')),length(file_nms)>1)
     file_nms(2:end)=file_nms(1);
    end
    for loop_nms=1:length(file_nms)
     file_nm=which(file_nms(loop_nms).name);
     H_vec(count_ht)=HT(1,loop_ht);
     T_vec(count_ht)=HT(2,loop_ht);
     amps_hld       =zeros(1,length(probes));
     eval(['load ' file_nm ' outputs;'])
     count_pb=1;
     for probe=probes
      out=outputs{count_pb+length(probes)};
      amps_hld(count_pb)=out.value;
      clear out
      count_pb=count_pb+1;
     end
     inc_amp(count_ht)=mean(amps_hld);
     inc_std(count_ht)=std(amps_hld);
     clear outputs count_pb amps_hld
     count_ht=count_ht+1;
    end
   end
   
   %%% Calculate ratio trans/inc
   
   trans_coeff     = trans_amp./inc_amp;
   if 0
    trans_coeff_err(1,:) = (trans_amp+trans_std)./(inc_amp-inc_std);
    trans_coeff_err(2,:) = (trans_amp-trans_std)./(inc_amp+inc_std);
   else
    trans_coeff_err(1,:) = (trans_amp+0.001)./(inc_amp-0.001);
    trans_coeff_err(2,:) = (trans_amp-0.001)./(inc_amp+0.001);
   end
   
   if EGY; trans_coeff=trans_coeff.^2; trans_coeff_err = trans_coeff_err.^2; end
   
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
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% HIGHER-ORDER HARMONIC PROPORTION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  elseif length(data_out.name)==2
   
   %% INCIDENT WAVE
   
   if strfind(t_meth,'inc')
    
    if COMM; cprintf(0.4*[1,1,1],['>>> incident waves\n']); end
    
    for loop_ht=ht_inds
     if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
       '; Hs=' num2str(HT(1,loop_ht)) ...
       ' ; eps=' num2str(HT(4,loop_ht))'\n']); end
     for probe=probes
      Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,Tpers,data_out,...
       loop_ht,probe,WaveType,'LHS',DO_FPLT,'none',DO_FDSP,file_pre);
      if DO_FPLT
       pause
       close all
      end
     end
    end
    
   elseif strfind(t_meth,'calib')
    
    if COMM; cprintf(0.4*[1,1,1],['>>> calibration waves\n']); end
    
    for loop_ht=ht_inds
     if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
       '; Hs=' num2str(HT(1,loop_ht)) ...
       ' ; eps=' num2str(HT(4,loop_ht)) '\n']); end
     for probe=probes
      Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),[],Tpers,data_out,...
       loop_ht,probe,WaveType,'RHS-calib',DO_FPLT,'none',DO_FDSP,file_pre);
      if DO_FPLT
       pause
       close all
      end
     end
    end
    
   end
   
   inc_i=zeros(length(dir([file_pre '*'])),1); inc_ii=inc_i;
   inc_i_std=inc_i; inc_ii_std=inc_ii;
   count_ht=1;
   for loop_ht=ht_inds
    file_nms=dir([file_pre sprintf('%03g',loop_ht) '*']);
    for loop_nms=1:length(file_nms)
     file_nm=which(file_nms(loop_nms).name);
     H_vec(count_ht)  =HT(1,loop_ht);
     T_vec(count_ht)  =HT(2,loop_ht);
     eps_vec(count_ht)=HT(4,loop_ht);
     amps_hld         =zeros(2,length(probes)); std_hld=amps_hld;
     eval(['load ' file_nm ' outputs;'])
     count_pb=1;
     for probe=probes
      out    ={outputs{count_pb}.value};
      out_std={outputs{count_pb}.std};
      amps_hld(1,count_pb)=out{1};
      amps_hld(2,count_pb)=out{2};
      std_hld(1,count_pb) =out_std{1};
      std_hld(2,count_pb) =out_std{2};
      clear out
      count_pb=count_pb+1;
     end
     inc_i(count_ht)       = mean(amps_hld(1,:));
     inc_ii(count_ht)      = mean(amps_hld(2,:));
     inc_i_std(count_ht)   = mean(std_hld(1,:)); %std(amps_hld(1,:)); %
     inc_ii_std(count_ht)  = mean(std_hld(2,:)); %std(amps_hld(2,:)); %
     clear outputs count_pb amps_hld std_hld
     count_ht=count_ht+1;
    end
   end
   
   inc_coeff = inc_ii./inc_i; clear inc_i inc_ii
   
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
   
   %% TRANSMITTED WAVE
   
   if COMM; cprintf(0.4*[1,1,1],'>>> transmitted waves\n'); end
   
   for loop_ht=ht_inds
    if COMM; cprintf(0.4*[1,1,1],['>>>> Tm=' num2str(HT(2,loop_ht)) ...
      '; Hs=' num2str(HT(1,loop_ht)) ...
      ' ; eps=' num2str(HT(4,loop_ht)) '\n']); end
    for probe=probes
     Main_AttnData(HT(2,loop_ht),HT(1,loop_ht),conc,Tpers,data_out,...
      loop_ht,probe,WaveType,'RHS',DO_FPLT,'none',DO_FDSP,file_pre,t_meth);
     if DO_FPLT
      pause
      close all
     end
    end
   end
   
   trans_i=zeros(length(dir([file_pre '*'])),1); trans_ii=trans_i;
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
      out={outputs{count_pb}.value};
      trans_i(count_ht) =trans_i(count_ht)+out{1};
      trans_ii(count_ht)=trans_ii(count_ht)+out{2};
      clear out
      count_pb=count_pb+1;
     end
     clear outputs count_pb
     count_ht=count_ht+1;
    end
   end
   trans_i = trans_i/length(probes); trans_ii = trans_ii/length(probes);
   
   trans_coeff = trans_ii./trans_i; clear trans_i trans_ii
   
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
   
  end
  
  %Generate Matrix and save it to Data/Gen folder with Conc and Type in
  %title
  T = T_vec;
  trans_std = trans_coeff_err;
  trans = trans_coeff.';
  H = H_vec;
  
  MatFile_NM = strcat('Data/Gen/Trans',int2str(conc),WaveType(1:3));
  save(MatFile_NM,'T','trans_std','trans','H');
  
  
  

return


