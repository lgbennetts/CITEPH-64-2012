% Main_WaveData
%
% generate some simple data to test WDM method
%
% INPUTS:
%
% [R,A]    = location of probes/staffs (A in rads)
% CREATEIT = create data
% RUNIT    = perform WDM on data (1=yes,0=no)
% run_num  = vector of indices of runs
% ns       = sampling frequency
% n        = 2^12
% f        = wave frequency in Hz
% th       = vector of wave angles wrt x-axis
% k0       = wavenumber
%
% VARIABLES:
%
% data     = array of wave elevations in time for each staff 
%
% written by L Bennetts Jan 2013 / Adelaide

function Main_WaveData

%% Prelims

path(path,'../../EXTRA_MATLAB_Fns'); % For colour

if ~exist('RUNIT','var'); RUNIT=1; end
if ~exist('CREATEIT','var'); CREATEIT=1; end
if ~exist('run_num','var'); run_num=1; end
if ~exist('TYP','var'); TYP=0; end

if CREATEIT
 
 if or(CREATEIT==1,CREATEIT==2)

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
 
 end
 
 % waves
 
 if ~exist('f','var'); f = 1; end
 if ~exist('k0','var'); k0 = (2*pi*f)^2/9.81; end

 for run=run_num
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if CREATEIT==1      % simple harmonic wave
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  th = 1*pi/6; 
  
  phi = exp(1i*k0*(cos(th)*X+sin(th)*Y));
 
  t=0;
  tm(1) = t;
  for loop_t = 1:3*n
   data(loop_t,:) = real(phi*exp(-2i*pi*f*t));
   t=t+(1/ns);
   tm(loop_t+1) = t;
  end
 
  description = ['simple harmonic wave: ang=' num2str(180*th/pi) ...
     '; f=' num2str(f) ' Hz; k=' num2str(k0) ' 1/m'];
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif CREATEIT==2     % superposition of harmonic waves
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  th = [-0,2]*pi/6; 
 
  phi = exp(1i*k0*(cos(th(1))*X+sin(th(1))*Y));
  for loop_a=2:length(th)
   phi = phi + exp(1i*k0*(cos(th(loop_a))*X+sin(th(loop_a))*Y));
  end
 
  t=0;
  tm(1) = t;
  for loop_t = 1:3*n
   data(loop_t,:) = real(phi*exp(-2i*pi*f*t));
   t=t+(1/ns);
   tm(loop_t+1) = t;
  end
 
  angs = num2str(180*th(1)/pi);
  for loop_a=2:length(th); angs = [angs ', ' num2str(180*th(loop_a)/pi)]; end
         
  description = [int2str(length(th)) ' simple harmonic waves: angs=' ...
     angs '; f=' num2str(f) ' Hz; k=' num2str(k0) ' 1/m'];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif CREATEIT==3 % data from expts LHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Tp=2; Hm=100;
  
  n = 2^12;
  
  dum_path = ['/Volumes/My_Passport/CITEPH_2012-053_CITEPH_WAVE_'...
    'PROPOGATION_IN_ICE-COVERED_SEAS/results_preliminary/'];
  
  c79_prams = conc79_testspecs();
  
  for loop=1:length(c79_prams)
   pers(loop)=10\c79_prams(loop).period;
   hts(loop) =10*c79_prams(loop).wave_height;
  end
  
  test = find(and(pers==Tp,hts==Hm));
  
  if isempty(test)
   cprintf('red',['Cant find Hm=' num2str(Hm) '; Tp=' num2str(Tp) '\n'])
   return
  end
 
  [xy_lhs,xy_rhs] = citeph_sensor_spots();
  dum_nms = dir([dum_path c79_prams(test).dirname '/houle_reg_*']);
  np = size(xy_lhs,1);
  count = length(dum_nms)-56+1; % 56 files(20 probes, 20 zoom, 16 accel) 
  
  X(1)=xy_lhs(1,1); Y(1)=xy_lhs(1,2);
  dum=load([dum_path c79_prams(test).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  for loop_xy=2:np
   X(loop_xy)=xy_lhs(loop_xy,1); Y(loop_xy)=xy_lhs(loop_xy,2);
   dum=load([dum_path c79_prams(test).dirname '/' dum_nms(count).name]);
   data(:,loop_xy) = dum(:,2)/100; clear dum 
   count=count+1;
  end
  ns = 1/tm(2);
  
  description = ['Oceanide expt: ' c79_prams(test).type ' waves;' ...
   ' wave maker side;' ...
   ' Hs=' num2str(10*c79_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c79_prams(test).period/10) ' [s]' ...
   ' f=' num2str(10/c79_prams(test).period) ' [Hz]'] ;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif CREATEIT==4 % data from expts RHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Tp=2; Hm=100;
  
  n = 2^12;
  
  dum_path = ['/Volumes/My_Passport/CITEPH_2012-053_CITEPH_WAVE_'...
    'PROPOGATION_IN_ICE-COVERED_SEAS/results_preliminary/'];
  
  c79_prams = conc79_testspecs();
  
  for loop=1:length(c79_prams)
   pers(loop)=10\c79_prams(loop).period;
   hts(loop) =10*c79_prams(loop).wave_height;
  end
  
  test = find(and(pers==Tp,hts==Hm));
  
  if isempty(test)
   cprintf('red',['Cant find Hm=' num2str(Hm) '; Tp=' num2str(Tp) '\n'])
   return
  end
 
  [xy_lhs,xy_rhs] = citeph_sensor_spots();
  dum_nms = dir([dum_path c79_prams(test).dirname '/houle_reg_*']);
  np = size(xy_lhs,1);
  count = length(dum_nms)-56+11; % 56 files(20 probes, 20 zoom, 16 accel) 
  
  X(1)=xy_rhs(1,1); Y(1)=xy_rhs(1,2);
  dum=load([dum_path c79_prams(test).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  for loop_xy=2:np
   X(loop_xy)=xy_lhs(loop_xy,1); Y(loop_xy)=xy_lhs(loop_xy,2);
   dum=load([dum_path c79_prams(test).dirname '/' dum_nms(count).name]);
   data(:,loop_xy) = dum(:,2)/100; clear dum 
   count=count+1;
  end
  ns = 1/tm(2);
  
  description = ['Oceanide expt: ' c79_prams(test).type ' waves;' ...
   ' beach side' ...
   ' Hs=' num2str(10*c79_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c79_prams(test).period/10) ' [s]' ...
   ' f=' num2str(10/c79_prams(test).period) ' [Hz]'] ;
  
 end
  
  %% Save data
    
  if run > 99
   eval(['save ../Fns/WDM/s13',int2str(run), ' X Y n ns data description tm TYP'])
  elseif run > 9
   eval(['save ../Fns/WDM/s130',int2str(run), ' X Y n ns data description tm TYP'])
  else
   eval(['save ../Fns/WDM/s1300',int2str(run), ' X Y n ns data description tm TYP'])
  end

 end % end run
end % end if CREATEIT
 
%% Run data
 
if RUNIT
 for run=run_num   
  WDM(run)
 end  
end
    
    
return   
    