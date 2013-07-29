function [time,data,file_list] = citeph_get_data(test_type,test_num,opt)
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

if strcmp(test_type,'calib');%%real or calibration
   fdir  = '/work/timill/CITEPH-data/calibration_waves/';
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
   fdir  = '/work/timill/CITEPH-data/results_preliminary/conc_79/';
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
   fdir  = '/work/timill/CITEPH-data/results_preliminary/conc_39/';
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


Nfiles   = length(dir([fdir,expt_dir '/*.dat']));

if strcmp(test_type,'calib');
   %% prev 20 = full eta series
   %% last 20 = zoom eta series
   nn    = [20 20];
   if strcmp(opt,'S')
      n0 = Nfiles+2 - sum(nn);
   else
      n0 = Nfiles+2 - nn(2);
   end
   nstep = 1;
else
   %% prev 20 = full eta series
   %% prev 20 = zoom eta series
   %% last 18 = [a_x,a_y,a_z] series
   nn = [20 20 16];
   if strcmp(opt,'S')
      n0    = Nfiles - sum(nn);
      nstep = 1;
   elseif strcmp(opt,'S_zoom')
      n0    = Nfiles - sum(nn(2:3));
      nstep = 1;
   elseif strcmp(opt,'Ax')
      n0    = Nfiles - nn(3);
      nvec  = [2 4 6 8 11 14];
   elseif strcmp(opt,'Ay')
      n0    = Nfiles - nn(3);
      nvec  = [1 5 7 9 12 15];
   else%% 'Az'
      %Nfiles,nn(3)
      n0    = Nfiles - nn(3);
      nvec  = [3 10 13 16];
   end
end

if strcmp(opt,'S')%%get full time series from all probes
   Nprobes     = 20;
   file_list   = cell(Nprobes,1);
   str0        = '000';
   for j=1:Nprobes
      jstr              = num2str(n0+j*nstep);
      nj                = length(jstr);
      lab               = str0;
      lab(end+1-nj:end) = jstr;
      %%
      fname          = [fdir,'/',expt_dir,'/',str1,str2,lab,'.dat']
      file_list{j}   = fname;
      %%
      A     = load(fname);
      if j==1
         time  = A(:,1);
         N     = length(time); 
         data  = zeros(N,Nprobes);
      end
      {data,A}
      data(:,j)   = A(:,2);
   end
elseif strcmp(opt,'S_zoom')%%zoomed timeseries
   Nprobes     = 20;
   file_list   = cell(Nprobes,1);
   str0        = '000';
   for j=1:Nprobes
      jstr              = num2str(n0+j*nstep);
      nj                = length(jstr);
      lab               = str0;
      lab(end+1-nj:end) = jstr;
      %%
      fname          = [fdir,'/',expt_dir,'/',str1,str2,lab,'.dat']
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
   str0        = '000';
   for j=1:Nprobes
      jstr              = num2str(n0+nvec(j));
      nj                = length(jstr);
      lab               = str0;
      lab(end+1-nj:end) = jstr;
      %%
      fname          = [fdir,'/',expt_dir,'/',str1,str2,lab,'.dat']
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
