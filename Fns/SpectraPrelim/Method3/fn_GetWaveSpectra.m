% function fn_JustFourierMovingWindow
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
% tm - time series times in seconds
% data - data at times (such as displacement in meters)
% Tpers - string containing function call, to be evaluated to give Tpers.
%
% General:
%
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function [Sn,Time,Periods,SnTimes] = fn_GetWaveSpectra(tm,disp,Tp,conc,PX,dt,smallesttimewindow,t0description,t1description,WaveType)

Tind = fn_Tind(conc,Tp/2,PX,WaveType); %Slower
% Tind = fn_Tind(conc,Tp,PX,WaveType);
%TindB = fn_Tind(conc,Tp,ProbeLocXY(1));
t0index = find(strcmp({Tind.description},t0description));
t0 = Tind(t0index).time;
t1index = find(strcmp({Tind.description},t1description));
t1 = Tind(t1index).time;  

MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
SI = MI - smallesttimewindow/2;
EI = MI + smallesttimewindow/2 -1;


% [timewindowlength,tm(SI),tm(EI)]

%% Fourier Transform
[Sn,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
Periods = 1./ff;

Time = [tm(SI), tm(EI)];


return





function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return


