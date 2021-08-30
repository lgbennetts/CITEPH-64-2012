function [xi,displacement,potential,R,T,f,ab,lambda,diff,D4Disp]=elastic_plate_modes(alpha,beta,gamma,h,L,n,N)
% [xi,displacement,potential,R,T]=elastic_plate_modes(alpha,beta,gamma,h,L,n,N)
% program to calculate the modes of a floating elastic plate
%--------------------------------------------------------------------------
% Variables:
% alpha  = frequency squared
% beta   = non-dimensional plate stiffness
% gamma  = non-dimensional plate mass
% h      =  water depth
% L      =  beam/plate half-length
% n      = number of modes
% N     =  number of calculation points along the plate is 2*N+1 (because
% we use Simpson's rule
%
% xi are the coefficients for the displacement, and displacement and
% potential and the values of the displacement and potential at the 2N+1 points
% evenly spaced between -L and L.
%
% xi = coefficients in the free modes 
% displacement = complex displacement
% potentiall = complex potential
% R = Reflection coefficient
% T = Transmission coefficient
% 
% Details can be found on http://www.wikiwaves.org/index.php/Green_Function_Methods_for_Floating_Elastic_Plates


% start of main

np = 2*N+1;                                                  


% Functions:
% beam_ev               ->  mode eigenvalues
% beam_em               ->  dry natural modes
% matrix_G_surface      ->  free-surface Green function
%--------------------------------------------------------------------------

% solve for the mode eigenvalues
lambda = beam_ev(n,L);

% define calculation points
x = -L:L/N:L;

% calculate dry natural modes
wn = beam_em(lambda,L,x);

% calculate free-surface Green function
G = matrix_G_surface(alpha,h,L,N);

% solve for the plate radiation potentials
phi_r = (eye(size(G))-alpha*G)\G*(-1i*sqrt(alpha))*wn.';

% calculate added mass and damping coefficients
% define simpsons coeff. matrix
dx = (2*L)/(np-1);
simp = 2*dx/3*ones(1,np);simp(2:2:np-1) = 4*dx/3;
simp(1) = dx/3;simp(end) = dx/3;
simp = diag(simp);

ab = transpose(wn*simp*phi_r);

k = dispersion_free_surface(alpha,0,h);
% calculate incident wave potential
phi_i =exp(-k*x);

% solve for the diffraction wave potential
phi_d = (eye(length(G))-alpha*G)\phi_i.';

% solve for the xi_n's
diff= wn*simp*phi_d; % calculate diffraction matrix
f = -1i*sqrt(alpha)*diff;

xi = (eye(n) + beta*diag(lambda.^4) - alpha*gamma*eye(n) + 1i*sqrt(alpha)*ab)\(-1i*sqrt(alpha)*diff); % calculated alpha_n

cond( 1/(-1i*sqrt(alpha))*eye(n) + (beta/(-1i)*sqrt(alpha) )*diag(lambda.^4) + sqrt(alpha)./1i*gamma*eye(n) - ab)


displacement = transpose(xi)*wn;
D4Disp= transpose(xi)*(diag(lambda.^4)*wn);
potential =  transpose(xi)*transpose(phi_r) + transpose(phi_d);



% we calculate the reflection and transmission coefficients.

R = -alpha/(tan(k*h) + k*h*sec(k*h)^2)*((exp(-k*x)).*(potential - 1i/sqrt(alpha)*displacement))*diag(simp);
T = 1 - alpha/(tan(k*h) + k*h*sec(k*h)^2)*((exp(k*x)).*(potential - 1i/sqrt(alpha)*displacement))*diag(simp);


% end of main

