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

function [Sn,mt,ff] = fn_JustFourierTransform(tm,data)


%% Fourier Transform

[Sn,mt,ff] = get_ft(tm,data);


return





%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%
%%%%%%%%%%%%%%%%%%


% function [Sn,ff,an,t_mean,xtra] = get_ft(time,displ)
%
% DESCRIPTION: fft in given time window
%              high pass filter by subtracting the mean;
%
% INPUTS:
%
% time      = sampling points in time
% displ     = time series of given quantity
% data_type = 1 (elevation) or 2 (acceleration - need to integrate F. series wrt time)
%             NO LONGER USED!!!!!!
%
% OUTPUTS:
%
% Sn     = spectrum corresponding to the positive frequencies ff (ff(1)=0);
% an     = magnitude of Fourier coefficients corresponding to 
%          [ff(1:end);-fmax;-flipud(ff(end:-1:2))] (ie no fftshift applied)
% ag     = angle corresponding to an
% t_mean = mean(time)
% xtra   = {Hs,Tp,T_m02}, T_m02=sqrt(m0/m2);
% 
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
%
%                N
% => displ(n) = sum a(fk)*exp(-2i*pi*(n-1)*(k-1)/N), 1 <= n <= N.
%               k=1
%
% where tn = time(1) + (n-1)*dt is discrete time  
% and   fk =           (k-1)*df is discrete frequency (df=1/N/dt)
%
% then
%
%                N
% a(fk) = (1/N) sum displ(tn)*exp(2i*pi*(n-1)*(k-1)/N), 1 <= k <= N.
%               n=1
%
% L Bennetts & T Williams 2013 / La Seyne & Adelaide

function [Sxx,t_mean,ff] = get_ft(time,displ)

displ = displ-mean(displ);

N        = length(time);

dt       = time(2)-time(1);                           %%time resolution
T        = time(N)+dt-time(1);                        %%record length
t_mean   = mean(time);
%%
df    = 1/T;                                          %%freq resolution
fs    = 1/dt;                                         %%sample rate
fmax  = fs/2;                                         %%Nyqvist freq (max possible freq)
%%
ff0   = (-N/2:N/2-1)'*df;
ff    = (0:N/2)'*df;
jj    = 1:N/2;
%%
an    = ifft(displ,N);
 
Sxx   = abs(an).^2;


Sneg  = [0;flipud(Sxx(jj+N/2))];                      %%-1,...,-N/2
Spos  = [Sxx(jj);0];                                  %%0,...N/2-1

Sxx   = Spos+Sneg;
 
ff(1) = []; Sxx(1) = [];

%an = sqrt(2*Sxx);

periods  = 1./ff;

%%

if 0
 tst_var  = [var(displ,1), sum(abs(an.^2)), sum(Sn)*df]
 pause
end

return

