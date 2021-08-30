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

function disp = fn_ExtractElevationFromEnergy(Sn,ff,tm)
% DETAILS:
%
% If displ(n) = eta(tn) for 1 <= n <= N
% &  an(k)    = a(fk)   for 1 <= k <= N,
%
% and we express
%
%            N
% eta(tn) = sum a(fk)*exp(-2i*pi*tn*fk), 1 <= n <= N.
%           k=1

%% Fourier Transform
% disp = zeros(1,length(tm));
% 
% for k= 1:length(ff)
%     disp = disp + (sqrt(2*Sn(k))*exp(-2i*pi*tm*ff(k)));
% end

disp = fft(sqrt(2*Sn),length(tm));


return


