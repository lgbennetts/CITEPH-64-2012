% Main_AttnData
%
% DESCRIPTION: either generate some simple data to test FFT/WDM analysis 
% methods or collect data from CITEPH experiments for analysis
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
%
% VARIABLES:
%
% data     = array of wave elevations in time for each staff 
%
% written by L Bennetts Jan 2013 / Adelaide

function Main_AttnData(Tp,Hm,conc,T_pers,data_out,run_num,probe,...
 CREATEIT,DO_PLOT,DO_SAVE,DO_DISP,file_pre)

if ~exist('Tp','var'); Tp=0.65; end
if ~exist('Hm','var'); Hm=20; end
if ~exist('conc','var'); conc=79; end

if ~exist('T_pers','var'); T_pers = 10; end

if ~exist('CREATEIT','var'); CREATEIT='LHS'; end

if ~exist('run_num','var');  run_num=1; end
if ~exist('RUNIT','var');    RUNIT='FFT'; end
if ~exist('TYP','var');      TYP=0; end
if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_PLOT','var');  DO_PLOT = 'Aspec'; end
if ~exist('DO_SAVE','var');  DO_SAVE = 0; end
if ~exist('DO_DISP','var');  DO_DISP = 1; end

if ~exist('probe','var'); 
 if strcmp(CREATEIT(2:3),'HS'); probe=10;
 elseif strcmp(CREATEIT(1:3),'ACC'); probe=1; end
end
if ~exist('zoom','var'); zoom=0; end

if ~exist('data_out','var') 
 data_out.name='amp-harmo-steady-1st';
 data_out.tint='t0=t0+4*Tp; t1=t0+10*Tp;'; 
end 

if ~exist('file_pre','var'); file_pre = 'Temp_data/a00'; end

dum_letters={'', 'a', 'b', 'c', 'd', 'e'};
  
if CREATEIT~=0
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if or(strcmp(CREATEIT,'plane1'),strcmp(CREATEIT,'plane2'))  % idealised
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % create an array of probes

 if 1
  A = [0 33:72:359]; dr=pi/180;  A = dr*A; 
  R = [0 0.25 0.25 0.25 0.25 0.25]; 
 else
  A = [0 120 240]; dr=pi/180;  A = dr*A; 
  R = [0.25 0.25 0.25];
 end
 
 X = R.*cos(A); Y = R.*sin(A);
  
 ns=4; 
 n=4096; 
 
 if ~exist('f','var'); f = 1; end
 if ~exist('k0','var'); k0 = (2*pi*f)^2/9.81; end
 
 Tp = 1/f; Hm = 1;
 
 data_type = 1;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif or(or(strcmp(CREATEIT,'LHS'),strcmp(CREATEIT,'RHS')),...
   strcmp(CREATEIT(1:3),'ACC'))                                % tests
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  basedir  = citeph_user_specifics;
  if isempty(basedir); return; end
  dum_path = [basedir 'attenuation_tests/Conc' int2str(conc) '/'];
  
  eval(['ex=exist(dum_path);'])
  
  if ex~=7
   cprintf('red','path does not exist\n')
   return
  end
  clear ex
  
  n = 2^12;
  
  eval(['c_prams = conc',int2str(conc),'_testspecs();'])
  
  pers = zeros(1,length(c_prams)); hts = pers;
  
  for loop=1:length(c_prams)
   if strcmp(c_prams(loop).type,'Regular')
    pers(loop)=10\c_prams(loop).period;
    hts(loop) =10*c_prams(loop).wave_height;
   end
  end
  
  test = find(and(pers==Tp,hts==Hm)); clear pers hts
  
  if isempty(test)
   cprintf('magenta',['Cant find Hm=' num2str(Hm) '; Tp=' num2str(Tp) '\n'])
   return
  end
  
  if length(test)>1
   cprintf('magenta',['Test repeated ' int2str(length(test)) ' times\n'])
   %t0=input('Which test(s)?');
   t0=1:length(test);
   test = test(t0); clear t0
   run_num=run_num+zeros(1,length(test));
  end
 
  [xy_lhs,xy_rhs] = citeph_sensor_spots();
  
  np = size(xy_lhs,1);
 
 end
 
 %%
 
 run_ct=1;
 
 for run=run_num
  
  %run=[int2str(run) dum_letters{run_ct}]; 
  
  dum_nms = ...
    dir([dum_path c_prams(test(run_ct)).dirname '/houle_reg_*']);
  
  if isempty(dum_nms)
   cprintf('magenta','Nothing in this folder\n')
   return
  end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if strcmp(CREATEIT,'plane1')      % simple harmonic wave
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  th = 0*pi/6; 
  
  t=0;
  for loop_t = 1:3*n
   data(loop_t,:) = cos(k0*(cos(th)*X+sin(th)*Y)-2*pi*f*t);
   t=t+(1/ns);
   tm(loop_t) = t;
  end
  
  data = (Hm/2)*data;
 
  description = ['simple harmonic wave: ang=' num2str(180*th/pi) ...
     '; f=' num2str(f) ' Hz; k=' num2str(k0) ' 1/m'];
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif strcmp(CREATEIT,'plane2')  % 2 simple harmonic waves
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  th = [-0,2]*pi/6; 
 
  phi = exp(1i*k0*(cos(th(1))*X+sin(th(1))*Y));
  for loop_a=2:length(th)
   phi = phi + exp(1i*k0*(cos(th(loop_a))*X+sin(th(loop_a))*Y));
  end
  phi = (Hm/2)*phi;
 
  t=0;
  for loop_t = 1:3*n
   data(loop_t,:) = real(phi*exp(-2i*pi*f*t));
   t=t+(1/ns);
   tm(loop_t) = t;
  end
 
  angs = num2str(180*th(1)/pi);
  for loop_a=2:length(th); angs = [angs ', ' num2str(180*th(loop_a)/pi)]; end
         
  description = [int2str(length(th)) ' simple harmonic waves: angs=' ...
     angs '; f=' num2str(f) ' Hz; k=' num2str(k0) ' 1/m'];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif strcmp(CREATEIT,'LHS') % data from expts LHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if conc==79
   count = length(dum_nms)-56+1+20*zoom; % 56 files(20 probes, 20 zoom, 16 accel) 
  elseif conc==39
   count = length(dum_nms)-54+1+20*zoom; % 54 files(20 probes, 20 zoom, 14 accel)
  end
  
  X(1)=xy_lhs(1,1); Y(1)=xy_lhs(1,2);
  dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  for loop_xy=2:np
   X(loop_xy)=xy_lhs(loop_xy,1); Y(loop_xy)=xy_lhs(loop_xy,2);
   dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
   data(:,loop_xy) = dum(:,2)/100; clear dum 
   count=count+1;
  end
  ns = 1/tm(2);
  
  X = X(probe); Y = Y(probe); data = data(:,probe);
  
  description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
   ' ' int2str(conc) '% conc;' ...
   ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  
  data_type = 1;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif strcmp(CREATEIT,'RHS') % data from expts RHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if conc==79
   count = length(dum_nms)-56+11+20*zoom; % 56 files(20 probes, 20 zoom, 16 accel) 
  elseif conc==39
   count = length(dum_nms)-54+11+20*zoom; % 54 files(20 probes, 20 zoom, 14 accel)
  end
  
  X(1)=xy_rhs(1,1); Y(1)=xy_rhs(1,2);
  dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  for loop_xy=2:np
   X(loop_xy)=xy_rhs(loop_xy,1); Y(loop_xy)=xy_rhs(loop_xy,2);
   dum=load([dum_path c_prams(test(run_ct)).dirname '/' dum_nms(count).name]);
   data(:,loop_xy) = dum(:,2)/100; clear dum 
   count=count+1;
  end
  ns = 1/tm(2);
  
  X = X(probe); Y = Y(probe); data = data(:,probe);
  
  description = ['Oceanide expt: ' c_prams(test(run_ct)).type ' waves;' ...
   ' ' int2str(conc) '% conc;' ...
   ' Hs=' num2str(10*c_prams(test(run_ct)).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test(run_ct)).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test(run_ct)).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  
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
  ' tm Tp TYP data_type data_out T_pers'])
 
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
  if     strcmp(RUNIT,'WDM')
   cprintf('green','>> WDM NEEDS WORK\n'); return
   WDM(run)
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
    