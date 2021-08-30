% function fn_InvestigateSpectra
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
%
% General:
%
% conc      = either 39 or 79 , otherwise calibration
% WaveType =  'Regular' or 'Irregular'
% TestName = name of test in fn_*_testspecs function for asociated
% experiment
% probe - probe to read, allow any from  1- 20 (10 on left, and 10 on right of probe)
%
%
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function [test_pram,test_WaveType,test_Success,calib_Success,calib_Name] ...
            = fn_FindTestAndCalib(conc,TestName)


%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Which user?
[user_number,user_name] = citeph_get_user();

if user_number == 0
    OTHER_USR=1;
end


% %% GENERAL
% %pick expriment
% if ~exist('conc','var') conc=79; end %39,79 or empty
% if ~exist('probe','var');     probe= 1 ;end; %probe can be 1,...,20
% if ~exist('TestName','var');     TestName= '6' ;end; 



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
% [xy_lhs,xy_rhs] = citeph_sensor_spots();
% xy_lhsrhs = [xy_lhs;xy_rhs];

if conc == 39 || conc == 79
    eval(['test_prams = conc',int2str(conc),'_testspecs();']);
else
    eval('test_prams = calib_testspecs();');
end    

%find a test matching conditions
IndexSatisfy =  find(ismember({test_prams.name},TestName)) ;

if isempty(IndexSatisfy)
    test_Success = 0;
    test_pram = [];
    test_WaveType = [];
    cprintf('red','No Experiment Matching Conditions \n');
    return
end
test_Success = 1;
test_pram = test_prams(IndexSatisfy);
test_WaveType = test_pram.type;
% ProbeLocXY,test_pram,test_WaveType,TestSuccess



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% Parameter Files / Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read in calibration
eval('c_prams = calib_testspecs();');
IndexSatisfyHs =  find([c_prams.wave_height]==test_pram.wave_height) ;


IndexSatisfyType =  find(strcmp({c_prams.type},test_pram.type));

IndexSatisfyTp =  find([c_prams.period]==test_pram.period) ;

IndexSatisfyHsTp = intersect(IndexSatisfyHs,IndexSatisfyTp);

IndexSatisfyHsTpType = intersect(IndexSatisfyHsTp,IndexSatisfyType);

if isempty(IndexSatisfyHsTpType)
    %Will have to use probes in front for calib
    calib_Success = 0;
    calib_Name = [];
    cprintf('red','No Calibration Matching Conditions');
    return
end
calib_Success = 1;
calib_Name = c_prams(IndexSatisfyHsTpType).name;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Generate  data matrices.
% 
% basedir  = citeph_user_specifics('final');
% 
% if conc == 39 || conc == 79
%     dum_path = [basedir '2-Wave_attenuation/Conc' int2str(conc) '/' WaveType '/'];
%     if strcmp(WaveType,'Regular')
%         dum_nms = ...
%         dir([dum_path c_prams(IndexSatisfy).dirname '/houle_reg_*']);
%     else
%         dum_nms = ...
%         dir([dum_path c_prams(IndexSatisfy).dirname '/houle_irr_*']);          
%     end
% else
%     dum_path = [basedir '1-Calibration/calibration_waves/' WaveType '/'];
%     if strcmp(WaveType,'Regular')
%         dum_nms = ...
%         dir([dum_path c_prams(IndexSatisfy).dirname '/calib_houle_reg_*']);
%     else
%         dum_nms = ...
%         dir([dum_path c_prams(IndexSatisfy).dirname '/calib_houle_irr_*']);          
%     end
% end
% 
% 
% 
% %Initial Sensor Index
% %conc - 39/79 or calibaration , regular
% if conc == 39
%     if strcmp(WaveType, 'Regular')
%         SensorInitialIndex = 13;
%     else
%         SensorInitialIndex = 7;
%     end
% elseif conc == 79
%     if strcmp(WaveType, 'Regular')
%         SensorInitialIndex = 13;
%     else
%         SensorInitialIndex = 7;
%     end 
% else
% %     if strcmp(WaveType, 'Regular')
% %         SensorInitialIndex = 13;
% %     else
% %         SensorInitialIndex = 7;
% %     end  
%     SensorInitialIndex = 7;
% end
% 
% ProbeLocXY = xy_lhsrhs(probe,:);
% 
% %Load appropriate data
% dum=load([dum_path c_prams(IndexSatisfy).dirname '/' ...
%     dum_nms(SensorInitialIndex - 1 +probe).name]);
% 
% tm = dum(:,1)/10;
% disp = dum(:,2)/100; clear dum;
% 
% Success = 1;

return
