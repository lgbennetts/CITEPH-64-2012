function x=runavf(x,fs,period);
%
% runavf.m;   1994-01-24 M. Donelan
%
% function x=runavf(x,fs,period);
%
%  This routine does a running average (time domain) on a time series
%  using the fft of a boxcar: sinx/x.
%  fs is sampling frequency in Hz.
%  period length of the boxcar, centred on the sample, in seconds.
%  
%
%

x=x(:);
lx=length(x);
if rem(lx,2)~=0;x=[x; x(lx)];end 
x=fft(x);
[n,m] = size(x);

	w=2*pi*(1:n/2)*fs/n;
	T=sin(w*period/2)./(w*period/2);

T=[1 T T(n/2-1:-1:1)]';
T=T*ones(1,m);
x=real(ifft(x.*T));
x(1,:) = x(2,:);
x(n,:) = x(n-1,:);
x=x(1:lx);
