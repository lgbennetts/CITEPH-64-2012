% Main_WaveData
%
% generate some simple data to test WDM method
%
% INPUTS:
%
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
% CREATEIT = create data (0 off, >0 identifier - see below)
% RUNIT    = perform analysis on data 
%            1=WDM 
%            2=MovingFFT
%            0=off
% TYP      = full directional analysis (1=yes, 0=no)
% DEL      = delete data after use (1=yes, 0=no)
%
% VARIABLES:
%
% data     = array of wave elevations in time for each staff 
%
% written by L Bennetts Jan 2013 / Adelaide

function Main_WaveData

%% Prelims

if ~exist('RUNIT','var'); RUNIT=2; end
if ~exist('CREATEIT','var'); CREATEIT=1; end
if ~exist('run_num','var'); run_num=1; end
if ~exist('TYP','var'); TYP=0; end
if ~exist('DEL','var'); DEL=1; end

if CREATEIT
 
 %%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if or(CREATEIT==1,CREATEIT==2)  % idealised
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
 
 Tp = 1/f; Hm = 2;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif or(CREATEIT==3,CREATEIT==4) % tests
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  basedir  = citeph_user_specifics;
  dum_path = [basedir '/results_preliminary/'];
  
  Tp=2; Hm=100;
  
  n = 2^12;
  
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
 
 end
 
 %%
 
 for run=run_num
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if CREATEIT==1      % simple harmonic wave
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif CREATEIT==2     % superposition of harmonic waves
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif CREATEIT==3 % data from expts LHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
   eval(['save ../Fns/WDM/s13',int2str(run), ' X Y n ns data description'...
    ' tm Tp TYP'])
  elseif run > 9
   eval(['save ../Fns/WDM/s130',int2str(run), ' X Y n ns data description'...
    ' tm Tp TYP'])
  else
   eval(['save ../Fns/WDM/s1300',int2str(run), ' X Y n ns data description'...
    ' tm Tp TYP'])
  end

 end % end run
end % end if CREATEIT
 
%% Run data
 
if     RUNIT==1
 for run=run_num   
  WDM(run)
 end  
elseif RUNIT==2
 for run=run_num   
  MovingFFT(run)
 end 
end % end if RUNIT

if DEL
 if run > 99
  eval(['delete ../Fns/WDM/s13',int2str(run) '.mat'])
 elseif run > 9
  eval(['delete ../Fns/WDM/s130',int2str(run) '.mat'])
 else
  eval(['delete ../Fns/WDM/s1300',int2str(run) '.mat'])
 end
end % end if DEL
    
    
return   
    