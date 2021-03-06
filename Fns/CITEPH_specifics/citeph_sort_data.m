function citeph_sort_data

tt = {'calib','c39','c79'};
for j=1:3
   test_type   = tt{j};
   citeph_sort_data0(test_type)
end

function citeph_sort_data0(test_type)

if nargin==0
   %test_type   = 'c39';
   %test_type   = 'c79';
   test_type   = 'calib';
end

if strcmp(test_type,'c79')
   Ntests   = 19;
elseif strcmp(test_type,'c39')
   Ntests   = 12;
elseif strcmp(test_type,'calib')
   Ntests   = 17;
end

%for test_num=2:16
for test_num=1:Ntests
   sort_data(test_type,test_num)
end

function sort_data(test_type,test_num)

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

basedir  = citeph_user_specifics;

if basedir(end)=='/'
 basedir(end)=[];
end

if strcmp(test_type,'calib');%%real or calibration
   fdir  = [basedir '/calibration/data_Wave_Calibration/calibration_waves/'];
   str1  = 'calib_houle_'; 
   %%
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_calib_prams(test_num);
   %%
   if strcmp(type,'Regular');%% reg/irreg waves
      str2  = 'reg_';
      nfil_vec = [20 20];
      str_vec  = {'wave_probes','wave_probes_zoom'};
   else
      str2  = 'irreg_';
      nfil_vec = [20 2];
      str_vec  = {'wave_probes','oceanide_spectra'};
   end
elseif strcmp(test_type,'c79')
   fdir  = [basedir '/results_preliminary/conc_79/'];
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
   fdir  = [basedir '/results_preliminary/conc_39/'];
   str1  = 'houle_'; 
   %%
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num);
   %%
   if strcmp(type,'Regular');%% reg/irreg waves
      str2     = 'reg_';
      nfil_vec = [20 20 14];
      %str_vec  = {'S','S_zoom','A'};
      str_vec  = {'wave_probes','wave_probes_zoom','accelerometers'};
   else
      str2     = 'irr_';
      nfil_vec = [20 14 2];
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
   dir2  = [fx_dir,'/',ftype];
   if ~exist(dir2)
      eval(['!mkdir ',dir2]);
   end
   Nmoved   = 0;
   while Nmoved<nfil_vec(j)
      nn    = nn-1;
      fname = DD(nn).name;
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
