% function fn_FindAndReadProbe
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% General: Find and extract probe information (time series - displacement)
%
% INPUTS:
%   Test information
%   conc - experiment
%   TestName - name of test in test_specs files
%   probe - number of probe to extract info from
%
% OUTPUTS:
%     ProbeLocXY - location of probe (x  [horizontal],y [vertical])
%     tm - time series
%     disp - displacement as function of tm
%     c_pram - parameters of the file, and test info
%     WaveType - type 'Regular' or 'Irregular'
%     Success - 1 if probe read successfully, otherwise 0
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by

function [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probe)


%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Which user?
[user_number,user_name] = citeph_get_user();

if user_number == 0
    OTHER_USR=1;
end


%% GENERAL
%pick expriment
if ~exist('conc','var') conc=79; end %39,79 or empty
if ~exist('probe','var');     probe= 1 ;end; %probe can be 1,...,20
if ~exist('TestName','var');     TestName= '6' ;end; 



%% DATA


% Consistency for other users:

if exist('OTHER_USR','var')
  cprintf('red','>>> USER NOT DETECTED, EITHER UPDATE INFO IN citeph_get_user and citeph_user_specifics to input data \n')  
  cprintf('red','    OR RELOAD TRANSMISSION MATRICES PROVIDED ON GITHUB \n') 
  error('NO USER DEFINED / LACK OF DATA REPOSITRY INFORMATION')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% Parameter Files / Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Read Param files 
[xy_lhs,xy_rhs] = citeph_sensor_spots();
xy_lhsrhs = [xy_lhs;xy_rhs];

if conc == 39 || conc == 79
    eval(['c_prams = conc',int2str(conc),'_testspecs();']);
else
    eval('c_prams = calib_testspecs();');
end    

%find a test matching conditions
IndexSatisfy =  find(ismember({c_prams.name},TestName)) ;

if isempty(IndexSatisfy)
    Success = 0;
    cprintf('red','No Experiment Matching Conditions \n');
    return
end

c_pram = c_prams(IndexSatisfy);
WaveType = c_pram.type;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Generate  data matrices.

basedir  = citeph_user_specifics('final');

if conc == 39 || conc == 79
    dum_path = [basedir '2-Wave_attenuation/Conc' int2str(conc) '/' WaveType '/'];
    if strcmp(WaveType,'Regular')
        dum_nms = ...
        dir([dum_path c_prams(IndexSatisfy).dirname '/houle_reg_*']);
    else
        dum_nms = ...
        dir([dum_path c_prams(IndexSatisfy).dirname '/houle_irr_*']);          
    end
else
    dum_path = [basedir '1-Calibration/calibration_waves/' WaveType '/'];
    if strcmp(WaveType,'Regular')
        dum_nms = ...
        dir([dum_path c_prams(IndexSatisfy).dirname '/calib_houle_reg_*']);
    else
        dum_nms = ...
        dir([dum_path c_prams(IndexSatisfy).dirname '/calib_houle_irr_*']);          
    end
end



%Initial Sensor Index
%conc - 39/79 or calibaration , regular
if conc == 39
    if strcmp(WaveType, 'Regular')
        SensorInitialIndex = 13;
    else
        SensorInitialIndex = 7;
    end
elseif conc == 79
    if strcmp(WaveType, 'Regular')
        SensorInitialIndex = 13;
    else
        SensorInitialIndex = 7;
    end 
else
%     if strcmp(WaveType, 'Regular')
%         SensorInitialIndex = 13;
%     else
%         SensorInitialIndex = 7;
%     end  
    SensorInitialIndex = 7;
end

ProbeLocXY = xy_lhsrhs(probe,:);

%Load appropriate data
dum=load([dum_path c_prams(IndexSatisfy).dirname '/' ...
    dum_nms(SensorInitialIndex - 1 +probe).name]);

tm = dum(:,1)/10;
disp = dum(:,2)/100; clear dum;

Success = 1;

return
