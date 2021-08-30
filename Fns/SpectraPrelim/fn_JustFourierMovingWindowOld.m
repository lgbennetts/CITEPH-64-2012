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

function [Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,data,Tp,T_window)


%% Fourier Transform

% get the length of the time windows - these are a bit arbitrary - but are
% from the paper:
% "Window width is chosen to balance time
% localization (narrow windows) and spectral resolution (wide windows, e.g. [36]). The target width
% in the example shown is ten target wave periods. (Target widths are three wave periods for the
% largest target wave periods.)
% ... (Target widths are three wave periods for the largest target wave
% periods.)"
% eval(Tpers);


N     = length(tm);
dt    = tm(2)-tm(1);
T     = tm(end)-tm(1);
fs    = 1/dt;
fmax  = fs/2;

%%

% T_window = Tpers*Tp;
df       = 1/T_window;
Nw       = 2^ceil(log2(T_window/dt));
T_window = dt*Nw;
t_window = (0:dt:(Nw-1))'*dt;

%%
Nf_S  = fmax*T_window;
%%
t_fac    = 0.99; % fraction of overlap between windows [1]
t_req    = T_window*(1-t_fac); % time step between windows [seconds]
n_shift  = floor(t_req/dt);
t_req    = n_shift*dt;
t_fac    = 1-t_req/T_window;
%%
Nwindows = floor((T-T_window)/t_req);
Sn_mat   = zeros(Nf_S,Nwindows);
an_mat   = zeros(Nf_S,Nwindows);

if 0
 %m0 = floor(100/T_window)
 m0 = 10
else
 m0 = 1;
end
t_int    = zeros(Nwindows,2);
tm_vec   = zeros(1,Nwindows);
% Tp_vec   = zeros(Nwindows);
% Hs_vec   = zeros(Nwindows);


for m=m0:Nwindows
 jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
 time0       = tm(jj);
 disp0       = data(jj);
 t_start     = time0(1);
 t_int(m,:)  = time0([1 end]);
 if m==m0 
  [Sn0,tm_vec(m),ff] = get_ft(time0,disp0);
 else
  [Sn0,tm_vec(m)]    = get_ft(time0,disp0);
 end
 Sn_mat(:,m) = Sn0; clear Sn0
end



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

%han window - apply window function and then multiply by factor
% HanWin = hanning(N);
% an    = ifft(HanWin.*displ,N);
% an    = sqrt(8/3)*an;
%does nothing really


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

