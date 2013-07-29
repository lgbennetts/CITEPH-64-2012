%% GODA_tankmodes.m
%% Author: Timothy Williams
%% Date:   20130718, 21:00:58 CEST

function [std_error,Tm,Rm,m_vec] = GODA_tankmodes_scan(period,xy,data,tank_width,GET_REFL)

if nargin==0
   DO_TEST           = 1;
   period            = 1.85;
   [xy_lhs,xy_rhs]   = citeph_sensor_spots();
   tank_width        = 16;%m
   %%
   xy       = xy_rhs;
   m_vec    = [0 3]';
   m_vec0   = m_vec;
   A        = [1 .55]';
   if 1
      B  = [];
   else
      B  = 0.3*A;
   end
   if isempty(B)
      GET_REFL = 0;
   else
      GET_REFL = 1;
   end
end

g  = 9.81;
om = 2*pi/period;
k  = om^2/g;
L  = tank_width;

M  = floor(k*L/pi)
if DO_TEST==1
   beta_m   = m_vec*pi/L;
   alp_m    = sqrt(k^2-beta_m.^2);
   data     = get_data(xy,A,B,k,alp_m,beta_m);
end

for m=1:M
   m_vec              = [0 m]';
   [err,Tm,Rm,mvec2] = GODA_tankmodes(period,xy,data,m_vec,tank_width,GET_REFL);
   std_error0(m)     = err(1);
end

[std_error,m]           = min(std_error0);
m_vec                   = [0 m]';
[std_error,T,R,m_vec2]  = GODA_tankmodes(period,xy,data,m_vec,tank_width,GET_REFL);

if DO_TEST==1;
   A,T
   B,R
   m_vec0,m_vec
end

function data  = get_data(xy,A,B,k,alp_m,beta_m)
%%
Nprobes     = size(xy,1);
Nm          = length(A);
data        = zeros(Nprobes,1);
if isempty(B)
   GET_REFL = 0;
else
   GET_REFL = 1;
end

if 1
   %% Do we need to add noise in time domain?
   %% I think think that Gaussian noise in the time domain
   %% corresponds to Gaussian noise in the frequency domain?
   noise_fac   = 0.05
   noise       = noise_fac*max(abs(A))*randn(Nprobes,1);
end

for i=1:Nprobes
   x  = xy(i,1);
   y  = xy(i,2);
   for j=1:Nm
      bet      = beta_m(j);
      alp      = alp_m(j);
      data(i)  = data(i)-k*A(j)*exp(1i*alp*x)*cos(bet*y);
      if GET_REFL==1
         data(i)  = data(i)-k*B(j)*exp(-1i*alp*x)*cos(bet*y);
      end
   end
end

data  = data+noise;
