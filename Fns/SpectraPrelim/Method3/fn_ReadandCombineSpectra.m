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

function [AvgSn,Period,SnTimes] = fn_ReadandCombineSpectra(conc,TestName,Tp,probes,t0description,t1description,TWindow)

if ~exist('conc','var') conc=1; end %39,79 or empty (1)
if ~exist('probes','var');     probes= 1:10 ;end; 
if ~exist('TestName','var');     TestName= '17' ;end; 


%SmallestTimeWindow - want to combine wave spectra, easiest way is just to
%make sure they're the same size
timewindowlengths = zeros(size(probes));
i = 1;


smallesttimewindow = min(timewindowlengths);
i = 1;
for probei = probes
 %Read Data in Calib 1, at probe
[ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probei);


[Sn,Time,Period] = fn_GetWaveSpectra(tm,disp,Tp,conc,ProbeLocXY(1),dt,smallesttimewindow,t0description,t1description,WaveType,TWindow);

SnMats(:,i) = Sn;
Times(:,i) = Time;
Periods(:,i) = Period;
i = i+1;
end

AvgSn = mean(SnMats,2);
Period = Periods(:,1);
SnTimes = [max(Times(1,:)),min(Times(2,:))];

% figure();
% plot(Periods,SnMats);
% hold on;
% plot(Periods(:,1),mean(SnMats,2),'--k');

return


function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

