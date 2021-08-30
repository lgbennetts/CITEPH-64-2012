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

function [TestNames,Success] = fn_FindTests(conc,Tp,Hs,WaveType)


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
if ~exist('conc','var') conc=39; end; %39,79 or empty
if ~exist('Tp','var');     Tp= 9.5 ;end; %probe can be 1,...,20
if ~exist('Hs','var');     Hs= NaN ;end;
if ~exist('WaveType','var');     WaveType= 'Regular' ;end;


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


if conc == 39 || conc == 79
    eval(['c_prams = conc',int2str(conc),'_testspecs();']);
else
    eval('c_prams = calib_testspecs();');
end    

%find a test matching conditions
if isnan(Tp)
    IndexSatisfyTp = 1: length(c_prams);
else
   IndexSatisfyTp =  find([c_prams.period] == 10*Tp) ; 
end

if isnan(Hs)
    IndexSatisfyHs = 1: length(c_prams);
else
   IndexSatisfyHs =  find([c_prams.wave_height] == Hs/10.0) ; 
end

if isnan(WaveType)
    IndexSatisfyWt = 1: length(c_prams);
else
   IndexSatisfyWt =  find(ismember({c_prams.type},WaveType)) ; 
end

% IndexSatisfyTp =  find([c_prams.period] == Tp) ;
% IndexSatisfyHs =  find([c_prams.wave_height] == Hs) ;
% IndexSatisfyWt =  find(ismember({c_prams.type},WaveType)) ;

IndexSatisfy = intersect(IndexSatisfyTp,IndexSatisfyHs);
IndexSatisfy = intersect(IndexSatisfy,IndexSatisfyWt);

if isempty(IndexSatisfy)
    Success = 0;
    cprintf('red','No Experiment Matching Conditions \n');
    return
end

ListNames = {c_prams.name};


TestNames = ListNames(IndexSatisfy);
Success = 1;

% IndexSatisfy =  find(ismember({c_prams.name},TestName)) ;





return
