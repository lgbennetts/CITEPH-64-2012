function Sample = SamplingDef(Signal)

% SampleParam = SamplingDef(SignalParam)
%
% This function defines the sampling properties of the input signal, 
% defined by the time-dependent function SignalT and the structure Signal. 
% This is designed to match the definition of the Matlab functions fft and 
% ifft.

%% FFT sampling

%%% Sampling interval
Sample.Dt=Signal.T0/10;

%%% Sampling frequency
Sample.fsamp=1/Sample.Dt;

%%% Signal length (in s)
Sample.T=Signal.t3;

%%% Number of time samples
Sample.Nt=floor(Sample.T/Sample.Dt);

%%% Time vector
Sample.t=(0:(Sample.Nt-1))*Sample.Dt;

%%% Number of frequency samples
Sample.Nsamp=ceil(Signal.t4/Sample.Dt);

%%% Time dependent input signal
Sample.Amps_t = Signal.SignalT(Sample.t,Signal);

%%% Sampling frequency interval
Sample.Df=Sample.fsamp/Sample.Nsamp;

%%% Frequency vector
Sample.f=(0:(Sample.Nsamp-1)/2)*Sample.Df;

%% Signal FFT

%%% Equivalent transform for positive frequencies
Sample.Amps_f    = 2*fft(Sample.Amps_t,Sample.Nsamp);
Sample.Amps_f(1) = Sample.Amps_f(1)/2;

%%% Shannon condition
Sample.Amps_f(length(Sample.f)+1:end) = [];

%% Input Signal Plots
% figure(3)
% subplot(2,1,1)
% plot(Sample.t,Sample.Amps_t)
% subplot(2,1,2)
% plot(2*pi*Sample.f,real(Sample.Amps_f))
% hold on 
% plot(2*pi*Sample.f,imag(Sample.Amps_f),'r')

%% More restricting filter used to neglect zeros frequencies

%%% Cut-off frequency
Sample.f_max = 1/Signal.T0+0.5;

%%%% Number of frequency samples
Sample.Nfilt = floor(Sample.f_max/Sample.Df);

%%% Updated frequency vector
Sample.f=(0:(Sample.Nfilt-1))*Sample.Df;

%%% Filter
Sample.Amps_f(Sample.Nfilt+1:end) = [];