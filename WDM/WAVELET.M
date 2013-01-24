

function [W,a]=wavelet(x,omin,omax,nv,ns);

%function [W,a]=wavelet(x,omin,omax,nv,ns);
%
%  Octave min is related to the lowest frequency such as 2^omin
%  For example Octave min = -2 --> .25 Hz
%  Octave max for the Highest frequency 2^omax
%  nv == number of voices
%  ns == sampling frequency [Hz]
%
% Bertrand Chapron's Wavelet Transform code.

%    Syntax W complex coeff. = wavelet ( x signal a column vector)
%
%Even though the concept is pretty wide for scale analysis these
%frequencies will be related to the peak in the Fourier domain
%

n=length(x);
f=fft(x);f=f(1:n/2+1);f(1)=f(1)/2;

ntot=(abs(omax-omin)+1)*nv;
for i=1:ntot,F(:,i)=f;end
k=1;


for i=omin:omax,
      for j=0:nv-1,
      a(k)=2^(i+j/nv);
      w(:,k)=momo(1/a(k),n,ns);
      k=k+1;
      end
end


W=ifft(2*F.*w,n);


