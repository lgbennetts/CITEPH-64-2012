% function fn_JustFourierMovingWindow
%
% DESCRIPTION: Extract spectra value at target period.
%         each spectra at target period is extracted from a single probe (average over frequency bins and between begin time and end time)
%         the spectract at target period is averaged over multiple probes 
%
% INPUTS:
%     tm - time series times in seconds
%     data - data at times (such as displacement in meters)
%     T_window - length of time window
%     Hz - Frequency of samples
%     Tp - Period to extraxt spectra value at
%
% OUTPUTS:
%     Sn_mat - matrix of spectra (ff x tm_vec)
%     tm_vec - mean time of windows
%     ff - frequency spectra is over
%     Nw - length of time window
% Jordan Pitt - Adelaide - 2021

function [Sn_mat,tm_vec,ff,Nw] = fn_JustFourierMovingWindow(tm,data,T_window,Hz,Tp)

N     = length(tm);
Nw = T_window;

olp_fac = 0.5;


HanWindow = hann(Nw);

%want all frequencies above Tp/10.0
datafilt = highpass(data,Tp/4.0,Hz);

[datawindowedmatrix,tm_vec] = createRollingWindow(datafilt,tm,Nw ,olp_fac);

[Sn_mat,ff] = periodogram(datawindowedmatrix,HanWindow,[], Hz);


return


function [Dwin, tm_vec] = createRollingWindow(vec,tm,wl,olp)
    l = length(vec);
    shift = max(1,floor(wl*(1-olp)));
    NReps = floor((l-wl)/ shift) + 1;
    
    for i = 1: NReps
        tmMat(1:wl,i) = tm((1+(i-1)*shift): (1 + (i-1)*shift) + (wl-1));
        Dwin(1:wl,i) = vec((1+(i-1)*shift): (1 + (i-1)*shift) + (wl-1)) ;
    end
    
    tm_vec = mean(tmMat,1);
 
return




