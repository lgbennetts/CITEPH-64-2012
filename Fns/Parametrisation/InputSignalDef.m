function Signal = InputSignalDef(g,H)

% Signal = InputSignalDef(g,H)
%
% This function defines the time dependent input signal of the wavemaker
% motion and returns a structure containing all its properties.

%% Frequency forcing

%%% Type of frequency forcing 'freq' or 'wlength'
type_force = 'freq';

if strcmp(type_force,'freq')
    %%% Frequency (in Hz)
    Signal.f0 = 0.5;
    
    %%% Frequency parameter
    Signal.kappa = (2*pi*Signal.f0).^2/g;
    
    %%% Period (in s)
    Signal.T0 = 1/Signal.f0;
    
    %%% Wavelength (in m)
    Signal.lam0 = 2*pi/CalcRealRoot_PWC([H, Signal.kappa], ...
        'FS_DispRel_PWC', 'UppLimReal_FS_PWC', 1e-16); 
elseif strcmp(type_force,'wlength')
    %%% Wavelength (in m)
    Signal.lam0 = 1;
    
    %%% Frequency (in Hz)
    Signal.f0 = FindFreq_FS([inf,inf,inf,inf,g],2*pi/Signal.lam0, H)/2/pi;
    
    %%% Period (in s)
    Signal.T0 = 1/Signal.f0;
    
    %%% Frequency parameter
    Signal.kappa = (2*pi*Signal.f0).^2/g;
end

%% Amplitude forcing

%%% Wavemaker transfer function
WM_TF = TF_Wavemaker(Signal.lam0,H);

%%% Type of amplitude forcing 'wave' or 'WM'
type_amp = 'wave';

if strcmp(type_amp,'wave')
    %%% Wave steepness (in %)
    Signal.eps = 1;
    
    %%% Steady-wave amplitude (in m)
    Signal.Amp0 = Signal.eps*Signal.lam0/200;
    
    %%% WM paddle amplitude (in m)
    Signal.AmpWM = WM_TF*Signal.Amp0;

elseif strcmp(type_amp,'WM')
    %%% WM paddle amplitude (in m)
    Signal.AmpWM = 0.1;
    
    %%% Steady-wave amplitude (in m)
    Signal.Amp0 = Signal.AmpWM/WM_TF;
end

%% Windows and Ramps
Signal.t0=0;        % Wavemaker start (in s).
Signal.t1=2;        % Transition ramp/steady (in s).
Signal.t2=50;       % Start of the decaying ramp (in s).
Signal.t3=55;       % Transition ramp/relaxation (in s).
Signal.t4=200;      % End of the relaxation time (in s).       

%% Time-dependent input wavemaker signal
% Linear ramp
Signal.SignalT= @(t,Signal) Signal.AmpWM.*(((t-Signal.t0).*...
    heaviside(t-Signal.t0)-(t-Signal.t1).*heaviside(t-Signal.t1))/...
    (Signal.t1-Signal.t0)+((t-Signal.t3).*heaviside(t-Signal.t3)-...
    (t-Signal.t2).*heaviside(t-Signal.t2))/(Signal.t3-Signal.t2)).*...
    cos(2*pi/Signal.T0*t);
