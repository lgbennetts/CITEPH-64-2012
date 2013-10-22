function [time,data,file_list,sensor_names] = citeph_get_data(test_type,test_num,opt,scale)
%% citeph_get_data.m
%% Author: Timothy Williams
%% Date:   20130723, 12:21:29 CEST
%% CALL: [time,data] = citeph_get_data(test_type,test_num)
%% test_type:  'c79', 'c39', 'calib'
%% test_num:   1,2,...
%% opt:        'S', 'S_zoom', 'Ax', 'Ay', 'Az'
%% scale:      'true' (measurement/basin scale), 'full' (lengths increased by 100, times increased by 10)

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
basedir  = citeph_user_specifics;

if strcmp(test_type,'calib');%%real or calibration
   fdir  = [basedir '/calibration/data_Wave_Calibration/calibration_waves/'];
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
   nn = [20 20 16];
   if strcmp(opt,'S')
      sub_dir  = '/wave_probes/';
      for j=1:20
         sensor_names{j}   = ['S',num2str(j)];
      end
   elseif strcmp(opt,'S_zoom')
      sub_dir  = '/wave_probes_zoom/';
      nstep = 1;
      for j=1:20
         sensor_names{j}   = ['S',num2str(j)];
      end
   elseif strcmp(opt,'Ax')
      sub_dir  = '/accelerometers/';
      nvec  = [2 4 6 8 11 14];
      for j=1:6
         sensor_names{j}   = ['A',num2str(j),'x'];
      end
   elseif strcmp(opt,'Ay')
      sub_dir  = '/accelerometers/';
      nvec  = [1 5 7 9 12 15];
      for j=1:6
         sensor_names{j}   = ['A',num2str(j),'y'];
      end
   else%% 'Az'
      %Nfiles,nn(3)
      sub_dir  = '/accelerometers/';
      nvec  = [3 10 13 16];
      jz    = [1,4:6];%%sensors that had the z-acceleration;
      for j=1:4
         sensor_names{j}   = ['A',num2str(jz(j)),'z'];
      end
   end
end
fx_dir   = [fdir,expt_dir,sub_dir];
DD       = dir([fx_dir, '/*.dat']);
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
   data     = data/fac;
end
