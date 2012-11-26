function WM_TF=TF_Wavemaker(lambda,H)

% WM_Amp=TF_Wavemaker(omega,lambda,H,Param)
%
% This function contains the inverse transfer function of the wavemaker 
% which gives the amplitude of the wavemaker AmpWM in terms of the wave 
% amplitude Amp0 that we want to generate. This of course depends also on 
% the water depth H and the frequency parameter kappa. The wavemaker shape  
% function must be called in this routine. 

%%% Travelling wavenumber
k0 = 2*pi/lambda;

%%% Innerproduct between the 2 travelling vertical modes
A0 = innerproduct_vertical(k0,k0,H)*weight_0_PWC(H, k0)^2;

%%% Definition of the wavemaker geometry
z = linspace(-H,0,100).';  % Discretised z-axis
f = paddle_shape(z,H);  % Generation of the discretized geometry.

%%% Inner-product vertical mode/paddle shape
IP = Trapezoidal_rule(cosh(k0*(z+H))*weight_0_PWC(H, k0).*f,-H,0);

%%% Wavemaker transfer function
WM_TF=abs(A0/(IP*weight_0_PWC(H, k0)*sinh(k0*H)));


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