function citeph_sort_data(test_type)

for test_num=4:19
%for test_num=1:19
   sort_data(test_type,test_num)
end

function sort_data(test_type,test_num)
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
   test_num    = 19;
%  opt         = 'S';
%  opt         = 'S_zoom';
%  opt         = 'Ax';
%  opt         = 'Ay';
%  opt         = 'Az';
end

if strcmp(test_type,'calib');%%real or calibration
   fdir  = '/work/timill/CITEPH-data/calibration_waves/';
   str1  = 'calib_houle_'; 
   %%
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_calib_prams(test_num);
   %%
   if strcmp(type,'Regular');%% reg/irreg waves
      str2  = 'reg_';
      nfil_vec = [20 20];%%[S S_zoom]
      str_vec  = {'S','S_zoom'};
   else
      str2  = 'irreg_';
      nfil_vec = [20 20];%%[S S_zoom]
   end
elseif strcmp(test_type,'c79')
   fdir  = '/work/timill/CITEPH-data/results_preliminary/conc_79/';
   str1  = 'houle_'; 
   %%
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num);
   %%
   if strcmp(type,'Regular');%% reg/irreg waves
      str2     = 'reg_';
      nfil_vec = [20 20 16];
      %str_vec  = {'S','S_zoom','A'};
      str_vec  = {'wave_probes','wave_probes_zoom','accelerometers'};
   else
      str2     = 'irr_';
      nfil_vec = [20 16 2];
      %str_vec  = {'S','A','oceanide_spectra'};
      str_vec  = {'wave_probes','accelerometers','oceanide_spectra'};
   end
else%% 'c39'
   fdir  = '/work/timill/CITEPH-data/results_preliminary/conc_39/';
   str1  = 'houle_'; 
   %%
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num);
   %%
   if strcmp(type,'Regular');%% reg/irreg waves
      str2     = 'reg_';
      nfil_vec = [20 20 16];
      %str_vec  = {'S','S_zoom','A'};
      str_vec  = {'wave_probes','wave_probes_zoom','accelerometers'};
   else
      str2     = 'irr_';
      nfil_vec = [20 16 2];
      %str_vec  = {'S','A','oceanide_spectra'};
      str_vec  = {'wave_probes','accelerometers','oceanide_spectra'};
   end
end

fx_dir   = [fdir,expt_dir];
DD       = dir([fx_dir '/*.dat']);
Nfiles   = length(DD);

nn       = Nfiles+1;
for j=length(str_vec):-1:1
   ftype = str_vec{j}
   dir2  = [fx_dir,'/',ftype]
   if ~exist(dir2)
      eval(['!mkdir ',dir2]);
   end
   Nmoved   = 0;
   while Nmoved<nfil_vec(j)
      nn    = nn-1;
      fname = DD(nn).name
      %eval(['!ln -s ' fx_dir '/' fname ' ' dir2]);
      eval(['!mv ' fx_dir '/' fname ' ' dir2]);
      Nmoved   = Nmoved+1;
   end
end

%%rest of files are analysis from Oceanide:
ftype = 'oceanide_analysis';
dir2  = [fx_dir,'/',ftype]
if ~exist(dir2)
   eval(['!mkdir ',dir2]);
end
eval(['!mv ' fx_dir '/*.dat ' dir2])
%eval(['!rm -r ' fx_dir '/S'])
%eval(['!rm -r ' fx_dir '/A'])
%eval(['!rm -r ' fx_dir '/S_zoom'])
