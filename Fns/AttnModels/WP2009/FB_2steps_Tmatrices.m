%% FB_2steps_Tmatrices.m
%% Author: Timothy Williams
%% Date:   20140226, 10:21:16 CET
function [R,A,B,T] = FB_2steps_Tmatrices(TM0,TM1,ka,inc_waves)

%%Transition matrices for LH edge
%%TM0 = {R0p,T0p,R0m,T0m}
R0p   = TM0{1};%%N1xN1
T0p   = TM0{2};%%N2xN1
R0m   = TM0{3};%%N2xN2
T0m   = TM0{4};%%N1xN2

%%Transition matrices for RH edge:
%%TM1 = {R1p,T1p,R1m,T1m}
if isempty(TM1)
   R1p   = R0m;%%N2xN2
   T1p   = T0m;%%N3xN2
   R1m   = R0p;%%N3xN3
   T1m   = T0p;%%N2xN3
else
   R1p   = TM1{1};%%N2xN2
   T1p   = TM1{2};%%N3xN2
   R1m   = TM1{3};%%N3xN3
   T1m   = TM1{4};%%N2xN3
end
%%
N1 = size(R0p,1);%%N LHS roots
N2 = size(T0p,1);%%N middle roots
N3 = size(T1p,1);%%N RHS roots

%%Default incident waves:
Ip = eye(N1,1);
Im = zeros(N3,1);
if exist('inc_waves')
   if ~isempty(inc_waves)
      Ip = inc_waves{1};
      Im = inc_waves{2};
   end
end

A0 = T0p*Ip;
B0 = T1m*Im;
R0 = R0p*Ip;
T0 = R1m*Im;
%%
ex = exp(1i*ka);
if 1%%set M
   M  = 20;
   jj = 1:M;
else
   jj = find(abs(ex)>1e-10);
end
%%
Dt    = diag(ex(jj));
A0t   = A0(jj);
B0t   = B0(jj);
R1p_t = R1p(jj,jj);
R0m_t = R0m(jj,jj);
%%
Bt = (eye(M,M)-R1p_t*Dt*R0m_t*Dt)\(B0t+R1p_t*Dt*A0t);
A  = A0+R0m(:,jj)*Dt*Bt;
At = A(jj);
B  = B0+R1p(:,jj)*Dt*At;
%%
R  = R0+T0m(:,jj)*Dt*Bt;
T  = T0+T1p(:,jj)*Dt*At;
