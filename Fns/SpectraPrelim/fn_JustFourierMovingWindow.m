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

function [Sn_mat,tm_vec,ff,Nw] = fn_JustFourierMovingWindow(tm,data,T_window,Hz,Tp)

% if ~exist('tm','var') tm = 0:0.01:1; end
% if ~exist('data','var') data = 1:0.01:2; end
% if ~exist('T_window','var') T_window = 10; end
% if ~exist('Hz','var') Hz = 10; end

N     = length(tm);


%%
% Nw       = 2^floor(log2(T_window));
Nw = T_window;

olp_fac = 0.5;

% dataCorr = data;

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
%         meanvec = mean(vec((1+(i-1)*shift): (1 + (i-1)*shift) + (wl-1)));
        tmMat(1:wl,i) = tm((1+(i-1)*shift): (1 + (i-1)*shift) + (wl-1));
        Dwin(1:wl,i) = vec((1+(i-1)*shift): (1 + (i-1)*shift) + (wl-1)) ;
    end
    
    tm_vec = mean(tmMat,1);
 
return




