function [time,data,file_list] = citeph_get_data_v2(test_type,test_num,opt)
%% citeph_get_data.m
%% Author: Timothy Williams
%% Date:   20130723, 12:21:29 CEST
%% CALL: [time,data] = citeph_get_data(test_type,test_num)
%% these are full scale not basin scale

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
basedir  = citeph_user_specifics;

if strcmp(test_type,'calib');%%real or calibration
   fdir  = [basedir '/calibration_waves/'];
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
   fdir  = [basedir '/results_preliminary/conc_79/'];
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
   fdir  = [basedir '/results_preliminary/conc_39/'];
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

if strcmp(opt(1),'A')
   root_dir = [fdir,expt_dir '/accelerometers'];
   DD0   = dir([root_dir '/*.dat']);
   %%
   if strcmp(test_type,'c79')
      if strcmp(opt,'Ax')
         nvec  = [2 4 6 8 11 14];
      elseif strcmp(opt,'Ay')
         nvec  = [1 5 7 9 12 15];
      else%% 'Az'
         nvec  = [3 10 13 16];
      end
   elseif strcmp(test_type,'c39')
      if strcmp(opt,'Ax')
         nvec  = [2 4 6 9 12];
      elseif strcmp(opt,'Ay')
         nvec  = [1 5 7 10 13];
      else%% 'Az'
         nvec  = [3 8 11 14];
      end
   end
   DD = DD0(nvec);
elseif strcmp(opt,'S')
   root_dir = [fdir,expt_dir '/wave_probes'];
   DD       = dir([root_dir '/*.dat']);
elseif strcmp(opt,'S_zoom')
   root_dir = [fdir,expt_dir '/wave_probes_zoom'];
   DD       = dir([root_dir '/*.dat']);
end
Nprobes  = length(DD);

if ~strcmp(opt,'S_zoom')%%get full time series from all probes
   file_list   = cell(Nprobes,1);
   for j=1:Nprobes
      fname          = [root_dir '/' DD(j).name];
      file_list{j}   = fname;
      %%
      A     = load(fname);
      if j==1
         time  = A(:,1);
         N     = length(time); 
         data  = zeros(N,Nprobes);
      end
      data(:,j)   = A(:,2);
   end
else%% 'S_zoom'
   file_list   = cell(Nprobes,1);
   for j=1:Nprobes
      fname          = DD(j).name;
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
end
