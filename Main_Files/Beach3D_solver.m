function S = Beach3D_solver(alpha, M, P, N, Rb, omega, xb, xi)

% Scat_b=Beach3D_solver(alpha,M,P,N,Rb,omega,xb,xi)
%
% S is the scattering matrix of the beach. It relates the incident wave to
% the reflected wave by a predefined reflexion coefficient. If cp is the
% vector of amplitudes of the incident wave and cm the vector of amplitude
% of the reflected wave. The condition can expressed as:
%
%                              cm = S*cp
%
% alpha is an (M+1)*(N+1)-length vector containing all the longitudinal 
% wavenumbers (propagating + evanescent).
% M+1 is the number of transversal modes.
% P is the number of propagating horizontal modes.
% N+1 is the number of vertical modes.
% omega is the radian frequency.
% xb is the location of the beach.
% xi is the reference of the wave incident on the beach. If it is not
% specified, we suppose the wave is generated from negative infinity and
% the reference is taken at the location of the beach.

if nargin==6
    xi=xb;
end

% Reflection coefficient
Rb = beach_TF(Rb, omega);

%% Note
% Each vector and matrix involved has length or size (N+1)*(M+1). The
% entries are organised so that the unknown amplitude vectors with entries
% a_nm are: a = [a_00 ... a_0M a_10 ... a1M a20 ... ... aNM].', where m and
% n are summation indices for the transversal and vertical modes,
% respectively.

%% Scattering matrix
% We consider a beach that only generates a travelling reflected one.

%%% Reflection matrix
R = zeros((M+1)*(N+1),1);
for m = 1:P
    R((m-1)*(N+1)+1) = Rb;
end

S=diag(R)*diag(exp(1i*alpha*(xb-xi)));