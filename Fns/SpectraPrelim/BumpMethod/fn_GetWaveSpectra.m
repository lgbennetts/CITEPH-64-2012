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

function [Sn,Time,Periods,SnTimes] = fn_GetWaveSpectra(tm,disp,Tp,conc,PX,TWindow,t0description,t1description,WaveType)

Tind = fn_Tind(conc,Tp/2,PX,WaveType); %Slower
% Tind = fn_Tind(conc,Tp,PX,WaveType);
%TindB = fn_Tind(conc,Tp,ProbeLocXY(1));
t0index = find(strcmp({Tind.description},t0description));
t0 = Tind(t0index).time;
t1index = find(strcmp({Tind.description},t1description));
t1 = Tind(t1index).time;  

[Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,TWindow);
Periods = 1./ff;

[~,jj0]=min(abs(tm_vec-t0));
[~,jj1]=min(abs(tm_vec-t1));

Sn = mean(Sn_mat(:,jj0:jj1), 2);
Time = [tm_vec(jj0), tm_vec(jj1)];



return





function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return


