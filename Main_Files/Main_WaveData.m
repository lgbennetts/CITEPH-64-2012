% Main_WaveData
%
% either generate some simple data to test FFT/WDM analysis methods
% or collect data from CITEPH experiments for analysis
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
%
% VARIABLES:
%
% data     = array of wave elevations in time for each staff 
%
% written by L Bennetts Jan 2013 / Adelaide

function Main_WaveData(RUNIT,CREATEIT,Tp,Hm,conc,probe,zoom)

%% Prelims

if ~exist('RUNIT','var'); RUNIT=2; end
if ~exist('CREATEIT','var'); CREATEIT='lhs'; end
if ~exist('run_num','var'); run_num=1; end
if ~exist('TYP','var'); TYP=0; end
if ~exist('DEL','var'); DEL=1; end

if ~exist('Tp','var'); Tp=2; end 
if ~exist('Hm','var'); Hm=100; end
if ~exist('conc','var'); conc=79; end
if ~exist('probe','var'); 
 if strcmp(CREATEIT(2:3),'HS'); probe=10;
 elseif strcmp(CREATEIT(1:3),'ACC'); probe=1; end
end
if ~exist('zoom','var'); zoom=0; end
  
if CREATEIT~=0
 
 %%
 
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
 
 Tp = 1/f; Hm = 2;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif or(or(strcmp(CREATEIT,'LHS'),strcmp(CREATEIT,'RHS')),...
   strcmp(CREATEIT(1:3),'ACC'))                                % tests
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  basedir  = citeph_user_specifics;
  dum_path = [basedir '/results_preliminary/'];
  
  eval(['ex=exist(dum_path);'])
  
  if ex~=7
   cprintf('red','path does not exist\n')
   return
  end
  clear ex
  
  n = 2^12;
  
  eval(['c_prams = conc',int2str(conc),'_testspecs();'])
  
  for loop=1:length(c_prams)
   pers(loop)=10\c_prams(loop).period;
   hts(loop) =10*c_prams(loop).wave_height;
  end
  
  test = find(and(pers==Tp,hts==Hm));
  
  if isempty(test)
   cprintf('red',['Cant find Hm=' num2str(Hm) '; Tp=' num2str(Tp) '\n'])
   return
  end
 
  [xy_lhs,xy_rhs] = citeph_sensor_spots();
  dum_nms = dir([dum_path c_prams(test).dirname '/houle_reg_*']);
  np = size(xy_lhs,1);
 
 end
 
 %%
 
 for run=run_num
 
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
  
  count = length(dum_nms)-56+1+20*zoom; % 56 files(20 probes, 20 zoom, 16 accel) 
  
  X(1)=xy_lhs(1,1); Y(1)=xy_lhs(1,2);
  dum=load([dum_path c_prams(test).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  for loop_xy=2:np
   X(loop_xy)=xy_lhs(loop_xy,1); Y(loop_xy)=xy_lhs(loop_xy,2);
   dum=load([dum_path c_prams(test).dirname '/' dum_nms(count).name]);
   data(:,loop_xy) = dum(:,2)/100; clear dum 
   count=count+1;
  end
  ns = 1/tm(2);
  
  X = X(probe); Y = Y(probe); data = data(:,probe);
  
  description = ['Oceanide expt: ' c_prams(test).type ' waves;' ...
   ' wave maker side;' ...
   ' Hs=' num2str(10*c_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  
  data_type = 1;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 elseif strcmp(CREATEIT,'RHS') % data from expts RHS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  dum_nms = dir([dum_path c_prams(test).dirname '/houle_reg_*']);
  np = size(xy_lhs,1);
  count = length(dum_nms)-56+11+20*zoom; % 56 files(20 probes, 20 zoom, 16 accel) 
  
  X(1)=xy_rhs(1,1); Y(1)=xy_rhs(1,2);
  dum=load([dum_path c_prams(test).dirname '/' dum_nms(count).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2)/100; clear dum
  count=count+1;
  for loop_xy=2:np
   X(loop_xy)=xy_rhs(loop_xy,1); Y(loop_xy)=xy_rhs(loop_xy,2);
   dum=load([dum_path c_prams(test).dirname '/' dum_nms(count).name]);
   data(:,loop_xy) = dum(:,2)/100; clear dum 
   count=count+1;
  end
  ns = 1/tm(2);
  
  X = X(probe); Y = Y(probe); data = data(:,probe);
  
  description = ['Oceanide expt: ' c_prams(test).type ' waves;' ...
   ' beach side' ...
   ' Hs=' num2str(10*c_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test).period) ' [Hz]; ' ...
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
  
  dum_nms = dir([dum_path c_prams(test).dirname '/houle_reg_*']);
  np = length(inds);
  count = length(dum_nms)-56+40; % 56 files(20 probes, 20 zoom, 16 accel) 
  
  cprintf('magenta','Need the accelerometer positions\n')
  
  X = 0; Y = 0;
  
  dum=load([dum_path c_prams(test).dirname '/' dum_nms(count+inds(1)).name]);
  tm = dum(:,1)/10; tm = tm(:);
  data(:,1) = dum(:,2); clear dum
  
  for loop_xy=2:np
   X(loop_xy)=0; Y(loop_xy)=0;
   dum=load([dum_path c_prams(test).dirname '/' dum_nms(count+inds(loop_xy)).name]);
   data(:,loop_xy) = dum(:,2); clear dum 
  end
  ns = 1/tm(2);
  
  X = X(probe); Y = Y(probe); data = data(:,probe);
  
  description = ['Oceanide expt: ' c_prams(test).type ' waves;' ...
   ' accelerometers' ...
   ' Hs=' num2str(10*c_prams(test).wave_height) ' [mm];' ...
   ' Tm=' num2str(c_prams(test).period/10) ' [s];' ...
   ' f=' num2str(10/c_prams(test).period) ' [Hz]; ' ...
   CREATEIT ' probe(s)=' int2str(probe)] ;
  
  data_type = 2;
  
 end % if CREATEIT=='string'
  
  %% Save data
    
  if run > 99
   eval(['save ../Fns/WDM/s13',int2str(run), ' X Y n ns data description'...
    ' tm Tp TYP data_type'])
  elseif run > 9
   eval(['save ../Fns/WDM/s130',int2str(run), ' X Y n ns data description'...
    ' tm Tp TYP data_type'])
  else
   eval(['save ../Fns/WDM/s1300',int2str(run), ' X Y n ns data description'...
    ' tm Tp TYP data_type'])
  end

 end % end run
end % end if CREATEIT
 
%% Run data
 
if     strcmp(RUNIT,'WDM')
 for run=run_num   
  WDM(run)
 end  
elseif strcmp(RUNIT,'FFT')
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
    