function [S, F] = Wavemaker3D_solver(X0, k, Z_weight, alpha, Y_weight, ...
    H, omega, M, N, Rw, W, xw, xi)

% [S,F]=Wavemaker2D_solver(X0,k,Z_weight,alpha,Y_weight,H,omega,M,N,Rw,W,xw,xi)
%
% S is the scattering matrix induced by the wavemaker and F is the wave
% generation vector, such that for a wave generated towards positive 
% infinity of amplitude ap and an incident wave coming from positive 
% infinity of amplitude am, we have:
%
%                            ap = S * am + F
%
% X0 is the amplitude of the wavemaker.
% k is the vector of wavenumbers (1 travelling + N evanescent).
% Z_weight is the vector of weighting coefficients of the vertical 
% modes in the open water domain (\cosh k_n(z+H)).
% H is the depth.
% alpha is an (M+1)*(N+1)-length vector containing all the longitudinal 
% wavenumbers (propagating + evanescent).
% Y_weight is the weighting coefficient of the transversal modes of the 
% wavetank (\cos \beta_m(z+H)).
% omega is the radian frequency of the paddle motion.
% M+1 is the number of waves propagating horizontally
% N+1 is the number of roots from the dispersion relation in open water 
% taken into account.
% Rw is the amplitude reflexion coefficient induced by the wavemaker.
% W is the width of the tank
% xw is the location of the paddle.
% xi is the horizontal reference for the wave incident on the wavemaker. 
% If it is not specified we suppose the wave is generated from positive 
% infinity and the reference is taken at the location of the paddle.


if nargin==12
    xi=xw;
end

%% Definition of the wavemaker geometry
N_w = 100;                  % Number of discrete values along the z-axis
z   = linspace(-H,0,N_w);   % z-coordinate vector
f   = paddle_shape(z,H);    % Generation of the discretized geometry.

%% Note
% Each vector and matrix involved has length or size (N+1)*(M+1). The
% entries are organised so that the unknown amplitude vectors with entries
% a_nm are: a = [a_00 ... a_N0 a_01 ... aN1 a02 ... ... aNM].', where m and
% n are summation indices for the transversal and vertical modes,
% respectively.

%% Definition of the input forcing vector F 
F=zeros((N+1)*(M+1),1);
for n=1:N+1
    F(n) = 1i*omega*X0*W*Y_weight*...
        Trapezoidal_rule(f.'.*cosh(k(n)*(z+H))*Z_weight(n),-H,0);
end

%% Evaluation of the scattering matrix S 
% The kinematic condition is expressed in integrated form as
%
%         IP_vert*IP_trans*(Phi_P_dx * ap + Phi_M_dx * am) = Fw.


%%% Matching in vertical direction (Free/Free In-Prod)
IP_vert=eye((N+1)*(M+1));
for n=1:N+1
    IP_vert(n,n) = ...
        innerproduct_vertical(k(n),k(n),H)*Z_weight(n)^2;
end

%%% Matching in transversal direction (Natural transversal modes In-Prod)
IP_trans=eye((N+1)*(M+1));
for n=1:N+1
    IP_trans(n,n) = W/2*Y_weight^2;
end


%%% Longitudinal eigenfunctions at wavemaker boundary
Phix_P_dx=diag(1i*alpha); 
Phix_M_dx=Rw*diag(-1i*alpha.*exp(-1i*alpha*(xi-xw))); 


%%% Scattering matrix
S=-Phix_P_dx\Phix_M_dx;

%%% We reexpress F to match the equation given in preambule
F=(IP_vert*IP_trans*Phix_P_dx)\F;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% sub-routines

function f = paddle_shape(z,H)

% f = paddle_shape(z,H)
%
% Shape function describing the geometry of a hinged-type wavemaker paddle.
% z defines the discretisation of the vertical domain
% H is the water depth

N_w=length(z);          
Zh=0.85*H;            % hinged point of the wavemaker localised at z=-Zh.
f=zeros(N_w,1);       % Initializing the paddle shape function.

for p=1:length(z)
    if z(p)<-Zh
        f(p)=0;
    else
        f(p)=(z(p)+Zh)/Zh;
    end
end

function out = innerproduct_vertical(k,kappa,h1)
% program to calculate 
% \int_{-h}^{-h+h1}\cosh \kappa\left( z+h\right) 
% \cosh k\left(z+h\right) dz

% first we test for the case when \kappa and k are equal (or nearly equal)
if abs((k-kappa)/h1) > 1e-6
    out = 1/2 * sinh((k - kappa)*h1) / (k - kappa) ...
        + 1/2 * sinh((k + kappa)*h1) / (k + kappa);
else
    if k == 0
        out = h1;
    else
        out = (1/2)*(sinh(2*k*h1)/(2*k) + h1);
    end
end

function out = Trapezoidal_rule(f,x_min,x_max)
% Trepezoidal_rule(f,x_min,x_max)
% 
% This routine calculates the integral of the sampled function f on the
% interval (x_min,x_max) using the trapezoidal rule.

if size(f,1)>1
    f=f.';
end

N = length(f);  % Number of samples
out = (x_max-x_min)/(2*(N-1))*sum([f(1) 2*f(2:N-1) f(end)]);