% Main_AttnData
%
% DESCRIPTION: collect data from CITEPH attenuation experiments for analysis
%              or generate some simple data to test FFT/WDM methods 
%
% INPUTS:
%
% conc  = concentration
% Hm    = wave height
% Tp    = period
% probe = select the probe (for Oceanide test data only) 
%         1-10 with 10 being the centre
%         use fn_whichprobe(side,probe,fig) to identify probe
% zoom  = used zoomed in data (0=off,1=on)
%
% VARIABLES:
% [R,A]    = location of probes/staffs (A in rads)
% run_num  = vector of indices of runs
% ns       = sampling frequency
% n        = 2^12
% f        = wave frequency in Hz
% th       = vector of wave angles wrt x-axis
% k0       = wavenumber
%
% FLAGS:
% 
% CREATEIT = create data (0 off, ~=0 identifier - see below)
% RUNIT    = perform analysis on data 
%            'WDM' 
%            'FFT' = moving FFT
%            0=off
% TYP      = full directional analysis (1=yes, 0=no)
% DEL      = delete data after use (1=yes, 0=no)
% data_type = 1 (displacement), 2 (acceleration)
% what_tests = 'prelim' or 'final' 
%
% VARIABLES:
%
% data     = array of wave elevations in time for each staff 
%
% written by L Bennetts Jan 2013 / Adelaide

function Main_AttnData(Tp,Hm,conc,Tpers,data_out,run_num,probe,WaveType,...
 CREATEIT,DO_PLOT,DO_SAVE,DO_DISP,file_pre,reps,RUNIT)

if ~exist('conc','var'); conc=79; end
if ~exist('wbin','var'); wbin=3; end
if ~exist('WaveType','var');     WaveType='Regular'; end

if or(~exist('Tp','var'),~exist('Tp','var')); 
 HT=fn_WhatTestData(conc,WaveType,0); HT=HT(:,3); end

if ~exist('Tp','var'); Tp=HT(2); end
if ~exist('Hm','var'); Hm=HT(1); end

if ~exist('Tpers','var'); Tpers='Tpers=fn_Tpers(Tp);'; end

if ~exist('CREATEIT','var'); CREATEIT='RHS'; end

if ~exist('run_num','var');  run_num=1; end
if ~exist('RUNIT','var');    RUNIT='FFT'; end 
if ~exist('TYP','var');      TYP=0; end
if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_PLOT','var');  DO_PLOT = 'Aspec-signal'; end
if ~exist('DO_SAVE','var');  DO_SAVE = 0; end
if ~exist('DO_DISP','var');  DO_DISP = 1; end
if ~exist('reps','var');     reps = 'all'; end %'calib'; end %

if ~exist('probe','var'); 
 if strcmp(CREATEIT(2:3),'HS'); probe=1;
 elseif strcmp(CREATEIT(1:3),'ACC'); probe=10; end
end
if ~exist('zoom','var'); zoom=0; end

if ~exist('data_out','var') 
 data_out.name={'amp-harmo-steady-1st'};
 data_out.tint='[t0,t1] = fn_tint(Tp,Tpers,t0,tvec);'; 
end 

if ~exist('file_pre','var');   file_pre = 'Temp_data/eg00'; end
if ~exist('what_tests','var'); what_tests = 'final'; end

dum_letters={'', 'a', 'b', 'c', 'd', 'e'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET UP DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if CREATEIT~=0
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if or(~isempty(strfind(CREATEIT,'plane1')),...
   ~isempty(strfind(CREATEIT,'plane2')))  % idealised
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [xy_lhs,xy_rhs] = citeph_sensor_spots();
  
 ns=2^6; 
 n =100*Tp*ns; %4096; 
 
 %Tp = 1;
 
 Param = ParamDef_Oceanide;

 fortyp='freq';

 Forcing = Force_def(Param.g, Param.bed, fortyp, 1/Tp);

 k0 = 2*pi/Forcing.lam0;
 
 Hm = Hm/100;
 
 data_type = 1;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif or(or(strcmp(CREATEIT,'LHS'),strcmp(CREATEIT,'RHS')),...
   strcmp(CREATEIT(1:3),'ACC'))                                % tests
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  basedir  = citeph_user_specifics(what_tests);
  if isempty(basedir); return; end
  if strcmp(what_tests,'prelim')
   dum_path = [basedir 'attenuation_tests/Conc' int2str(conc) '/'];
  elseif strcmp(what_tests,'final')
   dum_path = [basedir '2-Wave_attenuation/Conc' int2str(conc) '/' WaveType '/'];
  end
  
  eval(['ex=exist(dum_path);'])
  
  if ex~=7
   cprintf('red','>>> path does not exist\n')
   return
  end
  clear ex
  
  n = 2^12;
  
  eval(['c_prams = conc',int2str(conc),'_testspecs();'])
  
  pers = zeros(1,length(c_prams)); hts = pers;
  
  for loop=1:length(c_prams)
   if strcmp(c_prams(loop).type,WaveType)
    pers(loop)=10\c_prams(loop).period;
    hts(loop) =10*c_prams(loop).wave_height;
   end
  end
  
  test = find(and(pers==Tp,hts==Hm)); clear pers hts
  
  if isempty(test)
   cprintf('magenta',['>>> Cant 1 find Hm=' num2str(Hm) '; Tp=' num2str(Tp) '\n'])
   return
  end
  
  if length(test)>1
   cprintf('magenta',['>>> Test repeated ' int2str(length(test)) ' times\n'])
   %t0=input('Which test(s)?');
   if strfind(reps,'calib')
    t0=1; cprintf('magenta',['   ... 1st test only chosen\n'])
   else
    t0=1:length(test); cprintf('magenta',['   ... all tests chosen\n'])
   end
   test = test(t0); clear t0 reps
   run_num=run_num+zeros(1,length(test));
  end
 
  [xy_lhs,xy_rhs] = citeph_sensor_spots();
  
  np = length(probe); %np = size(xy_lhs,1);
 
 elseif strfind(CREATEIT,'calib')
  
  basedir  = citeph_user_specifics(what_tests);
  if isempty(basedir); return; end
  if strcmp(what_tests,'prelim')
   dum_path = [basedir 'calibration/data_Wave_Calibration/calibration_waves/' WaveType '/'];
  elseif strcmp(what_tests,'final')
   dum_path = [basedir '1-Calibration/calibration_waves/' WaveType '/'];
  end
  
  eval(['ex=exist(dum_path);'])
  
  if ex~=7
   cprintf('red','path does not exist\n')
   return
  end
  clear ex
  
  n = 2^12;
  
  eval(['c_prams = calib_testspecs();'])
  
  pers = zeros(1,length(c_prams)); hts = pers;
  
  for loop=1:length(c_prams)
   if strcmp(c_prams(loop).type,WaveType)
    pers(loop)=10\c_prams(loop).period;
    hts(loop) =10*c_prams(loop).wave_height;
   end
  end
  
  test = find(and(pers==Tp,hts==Hm)); clear pers hts
  
  if isempty(test)
   cprintf('magenta',['>>> Cant 2 find Hm=' num2str(Hm) '; Tp=' num2str(Tp) ' for calibration tests...\n'])
   cprintf('magenta',['   ... switching to LHS... \n'])
   Main_AttnData(Tp,Hm,conc,Tpers,data_out,run_num,probe,...
    'LHS',DO_PLOT,DO_SAVE,DO_DISP,file_pre,reps)
   return
  end
  
  conc     = [];
 
  [xy_lhs,xy_rhs] = citeph_sensor_spots();
  
%   np = size(xy_lhs,1);
  np = length(probe);  

 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 run_ct=1;
 
 for run=run_num
  
  %run=[int2str(run) dum_letters{run_ct}]; 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if strfind(CREATEIT,'plane1')     % simple harmonic wave
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strfind(CREATEIT,'LHS')
   X=xy_lhs(probe,1); Y(1)=xy_lhs(probe,2);
  elseif strfind(CREATEIT,'RHS')
   X=xy_rhs(probe,1); Y(1)=xy_rhs(probe,2);
  end
  
  th = 0*pi/6; 
  
  t=0;
  for loop_t = 1:n
   data(loop_t,:) = cos(k0*(cos(th)*X+sin(th)*Y)-2*pi*t/Tp);
   t=t+(1/ns);
   tm(loop_t) = t;
  end
  
  data = (Hm/2)*data;
 
  description = ['simple harmonic wave: ang=' num2str(180*th/pi) ...
     '; f=' num2str(1/Tp) ' Hz; k=' num2str(k0) ' 1/m'];
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif strfind(CREATEIT,'plane2')  % 2 simple harmonic waves
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if strfind(CREATEIT,'LHS')
   X=xy_lhs(probe,1); Y(1)=xy_lhs(probe,2);
  elseif strfind(CREATEIT,'RHS')
   X=xy_rhs(probe,1); Y(1)=xy_rhs(probe,2);
  end
  
  th = [0,1]*pi; 
  
  ab = (Hm/2)*[1,0.5];
 
  phi = ab(1)*exp(1i*k0*(cos(th(1))*X+sin(th(1))*Y));
  for loop_a=2:length(th)
   phi = phi + ab(loop_a)*exp(1i*k0*(cos(th(loop_a))*X+sin(th(loop_a))*Y));
  end
  
  t=0;
  for loop_t = 1:n
   data(loop_t,:) = real(phi*exp(-2i*pi*t/Tp));
   t=t+(1/ns);
   tm(loop_t) = t;
  end
 
  angs = num2str(180*th(1)/pi);
  for loop_a=2:length(th); angs = [angs ', ' num2str(180*th(loop_a)/pi)]; end
         
  description = [int2str(length(th)) ' simple harmonic waves: angs=' ...
     angs '; f=' num2str(1/Tp) ' Hz; k=' num2str(k0) ' 1/m'];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif strfind(CREATEIT,'LHS') % data from expts LHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isempty(conc)
      if strcmp(WaveType,'Regular')
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/calib_houle_reg_*']);
      else
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/calib_houle_irr_*']);          
      end
  else
      if strcmp(WaveType,'Regular')
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/houle_reg_*']);
      else
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/houle_irr_*']);          
      end
  end
  
  if isempty(dum_nms)
   cprintf('magenta','Nothing in this folder 1\n')
   return
  end
  
  if isempty(conc)
   description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
   ' calibration test;' ...
   ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  else
   description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
   ' ' int2str(conc) '% conc;' ...
   ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  end
  
  X=xy_lhs(probe,1); Y(1)=xy_lhs(probe,2);
 
  if isempty(conc) % calibration tests
   count = length(dum_nms)-40+0+20*zoom; % 40 files(20 probes, 20 zoom)
   if probe == 10; probe=19; end         % files in non-intuitive order
  elseif conc==79
   count = length(dum_nms)-56+0+20*zoom; % 56 files(20 probes, 20 zoom, 16 accel) 
  elseif conc==39
   count = length(dum_nms)-54+0+20*zoom; % 54 files(20 probes, 20 zoom, 14 accel 
  end
  
%   X(1)=xy_lhs(1,1); Y(1)=xy_lhs(1,2);
%   dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
%   tm = dum(:,1)/10; tm = tm(:);
%   data(:,1) = dum(:,2)/100; clear dum
%   count=count+1;
%   for loop_xy=2:np
%    X(loop_xy)=xy_lhs(loop_xy,1); Y(loop_xy)=xy_lhs(loop_xy,2);
%    dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
%    data(:,loop_xy) = dum(:,2)/100; clear dum 
%    count=count+1;
%   end
%   ns = 1/tm(2);
%   
%   X = X(probe); Y = Y(probe); data = data(:,probe);
%   
%   data_type = 1;

  for loop_xy=1:np
   dum=load([dum_path c_prams(test(run_ct)).dirname '/' ...
    dum_nms(count+probe(loop_xy)).name]);
   if loop_xy==1; tm = dum(:,1)/10; tm = tm(:); end
   data(:,loop_xy) = dum(:,2)/100; clear dum 
  end
  ns = 1/tm(2);
  
  data_type = 1;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif strfind(CREATEIT,'RHS') % data from expts RHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if isempty(conc)
      if strcmp(WaveType,'Regular')
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/calib_houle_reg_*']);
      else
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/calib_houle_irr_*']);          
      end
  else
      if strcmp(WaveType,'Regular')
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/houle_reg_*']);
      else
        dum_nms = ...
        dir([dum_path c_prams(test(run_ct)).dirname '/houle_irr_*']);          
      end
  end
  
  if isempty(dum_nms)
   cprintf('magenta','Nothing in this folder 2\n')
   return
  end
 
  if isempty(conc)
   description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
    ' calibration test;' ...
    ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
    ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
    ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
    CREATEIT ' probe(s)=' int2str(probe)] ;
  else
   description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
   ' ' int2str(conc) '% conc;' ...
   ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  end
  
  X=xy_rhs(probe,1); Y(1)=xy_rhs(probe,2);
  if strcmp(WaveType,'Regular')
      if isempty(conc) % calibration tests
       count = length(dum_nms)-40+9+20*zoom; % 40 files(20 probes, 20 zoom)
       if probe == 10; probe=11; end          % files in non-intuitive order
      elseif conc==79
       count = length(dum_nms)-56+10+20*zoom; % 56 files(20 probes, 20 zoom, 16 accel) 
      elseif conc==39
       count = length(dum_nms)-54+10+20*zoom; % 54 files(20 probes, 20 zoom, 14 accel)
      end
  else
      if isempty(conc) % calibration tests
       count = 7; % 40 files(20 probes, 20 zoom)
       if probe == 10; probe=11; end          % files in non-intuitive order
      elseif conc==79
       count = 7; % 56 files(20 probes, 20 zoom, 16 accel) 
      elseif conc==39
       count = 7; % 54 files(20 probes, 20 zoom, 14 accel)
      end
  end
  
%   X(1)=xy_rhs(1,1); Y(1)=xy_rhs(1,2);
%   dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
%   tm = dum(:,1)/10; tm = tm(:);
%   data(:,1) = dum(:,2)/100; clear dum
%   count=count+1;
%   for loop_xy=2:np
%    X(loop_xy)=xy_rhs(loop_xy,1); Y(loop_xy)=xy_rhs(loop_xy,2);
%    dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
%    data(:,loop_xy) = dum(:,2)/100; clear dum 
%    count=count+1;
%   end
%   ns = 1/tm(2);
%   
%   X = X(probe); Y = Y(probe); data = data(:,probe);

  for loop_xy=1:np
   dum=load([dum_path c_prams(test(run_ct)).dirname '/' ...
    dum_nms(count+probe(loop_xy)).name]);
   if loop_xy==1; tm = dum(:,1)/10; tm = tm(:); end
   data(:,loop_xy) = dum(:,2)/100; clear dum 
  end
  ns = 1/tm(2);
  
  data_type = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  elseif strcmp(CREATEIT(1:3),'ACC') % data from expts accels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strcmp(CREATEIT(4),'x')
   inds = [8,11,14,2,4,6];
  elseif strcmp(CREATEIT(4),'y')
   inds = [9,12,15,1,5,7];
  elseif strcmp(CREATEIT(4),'z')
   inds = [10,13,16,3];
  end
  
  dum_nms = dir([dum_path c_prams(test(run_ct)).dirname '/houle_reg_*']);
  np = length(inds);
  count = length(dum_nms)-56+40; % 56 files(20 probes, 20 zoom, 16 accel) 
  
  cprintf('magenta','Need the accelerometer positions\n')
  
  X = 0; Y = 0;
  
  dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count+inds(1)).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2); clear dum
  
  for loop_xy=2:np
   X(loop_xy)=0; Y(loop_xy)=0;
   dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count+inds(loop_xy)).name]);
   data(:,loop_xy) = dum(:,2); clear dum 
  end
  ns = 1/tm(2);
  
  X = X(probe); Y = Y(probe); data = data(:,probe);
  
  description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
   ' accelerometers' ...
   ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  
  data_type = 2;
  
 end % if CREATEIT=='string'
  
  %% Save data
    
  %file_nm=fn_get_filenm(run,file_pre);
  file_nm=[file_pre sprintf('%03g',run) dum_letters{run_ct}];
  
  eval(['save ', file_nm, '-in X Y n ns data description'...
  ' tm Tp TYP data_type data_out Tpers wbin probe'])
 
 run_ct=run_ct+1;
 
 clear data dum description

 end % end run
 
 clear run_ct
 
end % end if CREATEIT
 
%% Run data

if RUNIT
 run_ct=1;
 for run=run_num
  %run=[int2str(run) dum_letters{run_ct}];
  file_nm=[file_pre sprintf('%03g',run) dum_letters{run_ct}];
  if     strcmp(RUNIT,'WT')
   cprintf('green','>> WDM NEEDS WORK\n'); 
   WaveletTrans(file_nm,DO_PLOT,DO_SAVE,DO_DISP)
  elseif strcmp(RUNIT,'FFT')
   MovingFFT(file_nm,DO_PLOT,DO_SAVE,DO_DISP)
  end % end if RUNIT
  
  %if ~exist('file_nm','var'); file_nm=fn_get_filenm(run); end
  if DEL
   eval(['delete ' file_nm '-in.mat'])
  end
  run_ct=run_ct+1;
 end
 clear run_ct
end % END RUNIT
    
return   

%%% SUBFUNCTIONS %%%

function file_nm=fn_get_filenm(run,file_pre)

if ~exist('file_pre','var'); file_pre = 'Temp_data/a13'; end

if run > 99
  file_nm = [file_pre run];
 elseif run > 9
  file_nm = [file_pre '0' run];
 else
  file_nm = [file_pre '00' run];
end

return 
    