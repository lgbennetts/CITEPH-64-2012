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

function fn_SingleFourierTransform(conc, probe,TestName)


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
if ~exist('conc','var') conc=1; end %39,79 or empty (1)
if ~exist('probe','var');     probe= 1 ;end; 
if ~exist('TestName','var');     TestName= '14' ;end; 

if ~exist('wbin','var');   wbin =3; end


if ~exist('data_out','var')
 data_out.name={'amp-harmo-steady-1st'};
 data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
end

if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

%Read Data
[ProbeLocXY,tm,disp,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName, probe);

Tp = c_pram.period/10.0;
dt = tm(2) - tm(1);
[Energy,Peaks] = fn_ExtractAllPeaks(tm,disp,Tp,conc,ProbeLocXY(1),dt);

return


function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return

