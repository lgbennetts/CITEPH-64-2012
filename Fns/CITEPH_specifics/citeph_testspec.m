%% citeph_testspec.m
%% Author: Timothy Williams
%% Date:   20130715, 15:02:21 CEST

function [data,fs,freqs,wave_TH,Hs_ind] = citeph_testspec(num,xyz)

noise_fac   = 0.025

if 1%%choose from 39% tests
   T  = (6.5:3:18.5)';
   H  =  [2 3 4 4 4]';
else%%79% tests
   M  = [6.5   2.0
         8.0   2.0
         8.0   4.0
         9.5   3.0
         11.0  4.0
         12.5  4.0
         14.0  4.0
         14.0  8.0
         15.5  4.0
         17.0  4.0
         18.5  4.0
         20.0  4.0
         20.0  10.0];
   T  = M(:,1);
   H  = M(:,2);
end

%% 1:100 scale
fac      = 100;
jt       = 5;
wave_TH  = [T(jt) H(jt)];
%%
H     = H(jt)/fac;
om0   = 2*pi/T(jt);
g     = 9.81;
lam0  = 2*pi*g/om0^2;%%deep water
lam   = lam0/fac
om    = sqrt(2*pi*g/lam);
T     = 2*pi/om;
%%
wave_TH  = [wave_TH;T H]


T_length = 20*T;%%need ~20 cycles to get a wave
df       = 1/T_length;%% 
Nt       = 2^15;
dt       = T_length/Nt;
fs       = 1/dt;%%sample rate
f_max    = fs/2%% Nyqvist frequency
                %%  - smallest period is 2*dt, so max freq is fs/2

tt    = dt*(1:Nt)';
freqs = df*(-Nt/2:Nt/2-1)';%%will give negative frequencies too;
freqs = freqs(find(freqs>0));

%min(freqs)
%return;

if ~exist('num')
   num   = 1;
end

xx       = xyz(:,1);
yy       = xyz(:,2);
Nprobes  = length(xx);
ww       = zeros(Nt,Nprobes);
k        = 2*pi/lam;
switch num

case 1
   A     = H*[1/2 1/3 1/4]';
   ang   = pi/180*[0 30 60]';
   for n = 1:Nprobes
      x        = xx(n);
      y        = yy(n);
      for r=1:length(A)
         kx = k*cos(ang(r));
         ky = k*sin(ang(r));
         ww(:,n)  = ww(:,n)+...
            real( A(r)*exp(1i*(kx*x-om*tt))*sin(ky*y) );%%elevation
      end
   end
   noise = noise_fac*max(A)*randn(size(ww));

case 2%%sine wave
   A  = H/2;
   for n = 1:Nprobes
      x        = xx(n);
      y        = yy(n);
      ww(:,n)  = real( A*exp(1i*(k*x-om*tt)) );%%elevation
   end
   noise = noise_fac*A*randn(size(ww));
end
%%
data     = ww+noise;
Hs_ind   = 4*sqrt(var(data));
