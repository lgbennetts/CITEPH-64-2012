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
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function [CalibTestName,Factor,Tp,Hs,WaveType,Success]  = fn_FindCalibration(conc,TestName)

if ~exist('conc','var') conc =39; end %39,79 or empty (1)
if ~exist('TestName','var') TestName='7'; end 

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Which user?
[user_number,user_name] = citeph_get_user();

if user_number == 0
    OTHER_USR=1;
end

eval(['c_pramsconc = conc',int2str(conc),'_testspecs();']);

IndexSatisfy =  find(ismember({c_pramsconc.name},TestName)) ;
if isempty(IndexSatisfy)
    Success = 0;
    cprintf('red','No Experiment Matching Conditions \n');
    return
end

Tp = c_pramsconc(IndexSatisfy).period;
Hs = c_pramsconc(IndexSatisfy).wave_height;
WaveType = c_pramsconc(IndexSatisfy).type;

eval('c_prams = calib_testspecs();');
    

%find a test matching conditions
ITP =  find([c_prams.period]==Tp) ;
IType =  find(ismember({c_prams.type},WaveType)) ;
ITPTyp = intersect(ITP,IType);
if isempty(ITPTyp)
    Success = 0;
    cprintf('red','No Experiment Matching Conditions \n');
return
end
    
IHS =  find([c_prams.wave_height]==Hs) ;

ITPHS = intersect(ITPTyp,IHS);

if isempty(ITPHS)
    Factor = Hs/c_prams(ITPTyp).wave_height;
    Success = 1;
    CalibTestName     = c_prams(ITPTyp).name;
else
    Factor = 1;
    Success = 1;
    CalibTestName     = c_prams(ITPTyp).name;
end

Tp = Tp/10.0;
Hs = Hs/100.0;


return



