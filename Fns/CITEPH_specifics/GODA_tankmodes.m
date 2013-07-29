%% GODA_tankmodes.m
%% Author: Timothy Williams
%% Date:   20130718, 21:00:58 CEST

function [std_error,T,R,m_vec2] = GODA_tankmodes(period,xy,data,m_vec,tank_width,GET_REFL)

if nargin==0
   DO_TEST           = 1;
   period            = 1.85;
   [xy_lhs,xy_rhs]   = citeph_sensor_spots();
   tank_width        = 16;%m
   %%
   xy       = xy_rhs;
   m_vec0   = [0 3]';
   m_vec    = [0 3]';
   A        = [1 .55]';
   if 0
      B  = [];
   else
      B  = 0.3*A;
   end
   if isempty(B)
      GET_REFL = 0;
   else
      GET_REFL = 1;
   end
else
   DO_TEST  = 0;
end

g  = 9.81;
om = 2*pi/period;
k  = om^2/g;
L  = tank_width;

beta_m	 = m_vec*pi/L;
alp_m	 = sqrt(k^2-beta_m.^2);
if DO_TEST==1
   alp_m
   [data,std_error_act]  = get_data(xy,A,B,k,m_vec0,L);
end

%%discard m's that don't give propating modes
jim   = find(imag(alp_m)>0);
if isempty(jim)==0
   m_vec(jim)  = [];
   beta_m(jim) = [];
   alp_m(jim)  = [];
end

%%matrix: -rows correspond to diff probes (xy val's)
%%	  -col's correspond to diff m's
Nprobes	 = length(data);
Nm	 = length(m_vec);
MT       = zeros(Nprobes,Nm);
MR       = zeros(Nprobes,Nm);

for i=1:Nprobes
   x  = xy(i,1);
   y  = xy(i,2);
   for j=1:Nm
      bet      = beta_m(j);
      alp      = alp_m(j);
      MT(i,j)  = -k*exp(1i*alp*x)*cos(bet*y);
      MR(i,j)  = -k*exp(-1i*alp*x)*cos(bet*y);
   end
end

%%least squares fit with pinv (pseudo inverse)
if GET_REFL==0%%no reflections
   MR = [];
end

M  = [MT MR];
R  = pinv(M)*data;

if GET_REFL==0
   T        = R;
   R        = [];
   jx       = 2:Nm;
   [E,jj]   = sort(abs(T(jx)));
   jx2      = [1,jx(jj)];
   T        = T(jx2);
   m_vec2   = m_vec(jx2);
   M        = M(:,jx2);
else
   T        = R(1:Nm);
   R(1:Nm)  = [];
   jx       = 2:Nm;
   [E,jj]   = sort(abs(R(jx)).^2+abs(T(jx)).^2);
   jx2      = [1,jx(jj)];
   T        = T(jx2);
   R        = R(jx2);
   m_vec2   = m_vec(jx2);
   M        = [MT(:,jx2) MR(:,jx2)];
end

std_error   = zeros(Nm,1);
for m=1:Nm
   mm = 1:m;
   for j=1:Nprobes
      if GET_REFL==0
         std_error(Nm+1-m) = std_error(Nm+1-m)+...
            abs(data(j)-M(j,mm)*T(mm))^2/Nprobes;
      else
         std_error(Nm+1-m) = std_error(Nm+1-m)+...
            abs(data(j)-MT(j,mm)*T(mm)-MR(j,mm)*R(mm))^2/Nprobes;
      end
   end
end

if DO_TEST==1;
   m_vec2,m_vec
   A,T
   B,R
   std_error,std_error_act
end

function [data,std_error] = get_data(xy,A,B,k,m_vec,L)
%%
beta_m	 = m_vec*pi/L;
alp_m	 = sqrt(k^2-beta_m.^2);
%%
Nprobes  = size(xy,1);
Nm       = length(A);
data     = zeros(Nprobes,1);
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
   noise       = noise_fac*max(abs(A))*randn(Nprobes,2)*[1;1i];
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

data        = data+noise;
std_error   = sum(abs(noise).^2)/Nprobes;

