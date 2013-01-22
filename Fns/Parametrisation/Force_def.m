function Forcing = Force_def(g,H,type_freq,lam0)

% Forcing = Force_def()
%
% Setup of forcing parameters

%% Frequency forcing

%%% Type of frequency forcing 'freq' or 'wlength'
if ~exist('type_freq','var'); type_freq = 'freq'; end
if isempty(type_freq); type_freq = 'freqh'; end

if strcmp(type_freq,'freq')
    %%% Frequency (in Hz)
    Forcing.f = .5;
    
    %%% Frequency parameter
    Forcing.kappa = (2*pi*Forcing.f).^2/g;
    
    %%% Wavelength (in m)
    Forcing.lam0 = 2*pi/CalcRealRoot_PWC([H, Forcing.kappa], 'FS_DispRel_PWC', ...
        'UppLimReal_FS_PWC', 1e-16); 
elseif strcmp(type_freq,'wlength')
    %%% Wavelength (in m)
    if ~exist('lam0','var') 
     Forcing.lam0 = 1; 
    else
     if isempty(lam0) 
      Forcing.lam0 = 1; 
     else 
      Forcing.lam0 = lam0;
     end
    end
    
    %%% Frequency (in Hz)
    Forcing.f = FindFreq_FS([inf,inf,inf,inf,g],2*pi/Forcing.lam0, H)/2/pi;
    
    %%% Frequency parameter
    Forcing.kappa = (2*pi*Forcing.f).^2/g;
elseif strcmp(type_freq,'waveno')
    %%% Wavelength (in m)
    if ~exist('lam0','var') 
     Forcing.lam0 = 1; 
    else
     if isempty(lam0) 
      Forcing.lam0 = 1; 
     else 
      Forcing.lam0 = 2*pi/lam0;
     end
    end
    
    %%% Frequency (in Hz)
    Forcing.f = FindFreq_FS([inf,inf,inf,inf,g],2*pi/Forcing.lam0, H)/2/pi;
    
    %%% Frequency parameter
    Forcing.kappa = (2*pi*Forcing.f).^2/g;    
end

%% Amplitude forcing

%%% Wavemaker transfer function
WM_TF = TF_Wavemaker(Forcing.lam0,H);

%%% Type of amplitude forcing 'wave' or 'WM'
type_amp = 'wave';

if strcmp(type_amp,'wave')
    %%% Wave steepness (in %)
    Forcing.eps = 1;
    
    %%% Steady-wave amplitude (in m)
    Forcing.Amp0 = Forcing.eps*Forcing.lam0/200;
    
    %%% WM paddle amplitude (in m)
    Forcing.AmpWM = WM_TF*Forcing.Amp0;
elseif strcmp(type_amp,'WM')
    %%% WM paddle amplitude (in m)
    Forcing.AmpWM = 0.1;
    
    %%% Steady-wave amplitude (in m)
    Forcing.Amp0 = Forcing.AmpWM/WM_TF;
end