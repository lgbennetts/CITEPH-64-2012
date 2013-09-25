% function Main_SingleDiskData
% 
% INPUTS:
%
% Hm    = wave height
% Tp    = period
% probe = select the probe (for Oceanide test data only)
%         1-10 with 10 being the centre
%         use fn_whichprobe(side,probe,fig) to identify probe
%         or, for the rigid body modes,
%         'surge'  amplitude in x-dir
%         'sway'   amplitude in y-dir
%         'heave'  amplitude in z-dir
%         'roll'   angle about x-axis
%         'pitch'  angle about y-axis
%         'yaw'    angle about z-axis
%
% VARIABLES:
% [R,A]    = location of probes/staffs (A in rads)
% run_num  = vector of indices of runs
% ns       = sampling frequency
% n        = 2^12
% f        = wave frequency in Hz
%
% FLAGS:
%
% CREATEIT  = create data (0 off, ~=0 identifier - see below)
% RUNIT     = perform analysis on data
%             'WDM'
%             'FFT' = moving FFT
%             0=off
% TYP       = full directional analysis (1=yes, 0=no)
% DEL       = delete data after use (1=yes, 0=no)
% data_type = 1 (displacement), 2 (acceleration), 3 (degrees)
%
% VARIABLES:
%
% data     = array of wave elevations in time for each staff
%
% written by L Bennetts Aug 2013 / Adelaide

function Main_SingleDiskData(Tp,Hm,T_pers,data_out,...
 run_num,probe,CREATEIT,DO_PLOT,DO_SAVE,DO_DISP,file_pre)

if ~exist('Tp','var'); Tp=1.85; end
if ~exist('Hm','var'); Hm=40; end
if ~exist('T_pers','var'); T_pers = 10; end

if ~exist('RUNIT','var');    RUNIT='FFT'; end
if ~exist('CREATEIT','var'); CREATEIT='WP'; end
if ~exist('run_num','var');  run_num=1; end
if ~exist('TYP','var');      TYP=0; end
if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_PLOT','var');  DO_PLOT = 'Aspec-signal'; end
if ~exist('DO_SAVE','var');  DO_SAVE = 0; end
if ~exist('DO_DISP','var');  DO_DISP = 1; end

if ~exist('probe','var');
 if strcmp(CREATEIT,'WP'); probe=2;
 elseif strcmp(CREATEIT,'RBM'); probe='surge'; end
end

if ~exist('data_out','var') 
 data_out.name='amp-harmo-steady-1st';
 data_out.tint='t0=t0+4*Tp; t1=t0+10*Tp;'; 
end 

if ~exist('file_pre','var'); file_pre = 'Temp_data/s00'; end

if CREATEIT~=0

%% Find file

basedir  = citeph_user_specifics;
dum_path = [basedir 'single_floe/'];

eval(['ex=exist(dum_path);'])

if ex~=7
 cprintf('red','path does not exist\n')
 return
end
clear ex

n = 2^12;

eval('c_prams = singlefloe_testspecs();')

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
 t0=input('Which test? ','s');
 eval(['test = test(' t0 '); clear t0'])
end

dum_nms = dir([dum_path c_prams(test).dirname '/drag_waves_*.dat']);

if isempty(dum_nms)
 cprintf('magenta','Nothing in this folder\n')
 return
end

%%

for run=run_num
 
 if strcmp(CREATEIT,'WP')
  
  xy_lhs = citeph_sensor_spots();
  
  np = size(xy_lhs,1);
  
  count = 3+probe; % wave probes raw data in files 4-13
  
  X(1)=xy_lhs(1,1); Y(1)=xy_lhs(1,2);
  dum=load([dum_path c_prams(test).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  
  ns = 1/tm(2);
  
  description = ['Oceanide expt: ' ...
   ' single floe;' ...
   ' Hs=' num2str(10*c_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test).period) ' [Hz]; ' ...
   CREATEIT ' probe=' int2str(probe)] ;
  
  data_type = 1;
  
 elseif strcmp(CREATEIT,'RBM')
  
  count = 13; % Rigid body mode data in files 14-19
  
  if strcmp(probe,'surge')
   prb_str=probe; probe=1; data_type=1;
  elseif strcmp(probe,'sway')
   prb_str=probe; probe=2; data_type=1;
  elseif strcmp(probe,'heave')
   prb_str=probe; probe=3; data_type=1;
  elseif strcmp(probe,'roll')
   prb_str=probe; probe=4; data_type=3;
  elseif strcmp(probe,'pitch')
   prb_str=probe; probe=5; data_type=3;
  elseif strcmp(probe,'yaw')
   prb_str=probe; probe=6; data_type=3;
  else
   cprintf('green','>>> ''probe'' not recognised\n'); return
  end
  
  X=0; Y=0;
  dum=load([dum_path c_prams(test).dirname '/' dum_nms(count+probe).name]);
  tm = dum(:,1)/10; tm = tm(:);
  if probe<=3; data(:,1) = dum(:,2)/100; clear dum; 
  else data(:,1) = dum(:,2); clear dum; end
  
  ns = 1/tm(2);
  
  description = ['Oceanide expt: ' ...
   ' single floe;' ...
   ' Hs=' num2str(10*c_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test).period) ' [Hz]; ' ...
   CREATEIT ' ' prb_str];
  
 end % end if CREATEIT
 
 %% Save data
 
 file_nm=fn_get_filenm(run,file_pre);
 
 eval(['save ', file_nm, '-in X Y n ns data description'...
  ' tm Tp TYP data_type data_out T_pers'])

 end % end run
 
end % end if CREATEIT

%% Run data

if RUNIT~=0
 
 if     strcmp(RUNIT,'WDM')
  for run=run_num
   cprintf('green','>> WDM NEEDS WORK\n'); return
   WDM(run)
  end
 elseif strcmp(RUNIT,'FFT')
  for run=run_num
   MovingFFT(file_nm,DO_PLOT,DO_SAVE,DO_DISP)
  end
 end % end if RUNIT
 
 if ~exist('file_nm','var'); file_nm=fn_get_filenm(run); end
 if DEL
  eval(['delete ' file_nm '-in.mat'])
 end
 
end % end if RUNIT
 

return

%%% SUBFUNCTIONS %%%

function file_nm=fn_get_filenm(run,file_pre)

if ~exist('file_pre','var'); file_pre = 'Temp_data/s13'; end

if run > 99
  file_nm = [file_pre int2str(run)];
 elseif run > 9
  file_nm = [file_pre '0' int2str(run)];
 else
  file_nm = [file_pre '00' int2str(run)];
end

return 
