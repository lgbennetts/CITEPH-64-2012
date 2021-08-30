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

function fn_StatisticsOverTime()

if ~exist('conc','var') conc=39; end %39,79 or empty (1)
if ~exist('probeT','var');     probeT= 11 ;end; 
if ~exist('probeR','var');     probeR= 9 ;end; 
if ~exist('TestName','var');     TestName= '2b' ;end; 
if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

if ~exist('wbin','var');     wbin= 3 ;end; 

close all;

%Reflection Probe

[ProbeLocXY,tmR,dispR,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName, probeR);

Tp = c_pram.period / 10.0;
Hm = c_pram.wave_height /100.0;

Description = ['Concentration ',num2str(conc) '%  , ProbeR : ', num2str(probeR),' ProbeT : ', num2str(probeT),'  Tp =', num2str(Tp), '[Hz]  Hs =', num2str(Hm), '[m]'  ];


T_window = 100*Tp;

[Sn_matR,tm_vecR,ffR] = fn_JustFourierMovingWindow(tmR,dispR,Tp,T_window);
periodR = 1./ffR;


%Transmission Probe

[ProbeLocXY,tmT,dispT,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName, probeT);
[Sn_matT,tm_vecT,ffT] = fn_JustFourierMovingWindow(tmT,dispT,Tp,T_window);
periodT = 1./ffT;

%Raw Signal
figure();
plot(tmR,dispR,'-k','DisplayName','Reflection Probe')
hold on;
plot(tmT,dispT,'-r','DisplayName','Transmission Probe')
hold off;
title(['Raw Signal ' , Description]);
xlabel('time (s)')
ylabel('displacement (m)')



return



