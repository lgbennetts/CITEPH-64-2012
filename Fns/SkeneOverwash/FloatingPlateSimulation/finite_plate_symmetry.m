function [a_s,b_s,a_a,b_a] = finite_plate_symmetry(alpha,beta,gamma,h,N,L,theta,nu)
% [a_s,b_s,a_a,b_a] = finite_plate_symmetry(alpha,beta,gamma,h,N,theta)
% 
% this program solves the finite plate problem by an
% eigenfunction matching method. 
% alpha, beta and gamma are non-dimensionnal parameters depending on the
% plate (details can be found on the wikiwave site).
% h is the water depth, 2L is the plate length, N is the number of real roots taken when solving
% the dispersion equations, and theta is the wave incident angle.
%
% a is the vector of coefficients for the eigenfunctions
% cos(k_n(z+h)/cos(k_n h) where k_n are the roots of the dispersion
% equation for a free surface
%
% b is the vector of coefficients for the eigenfunctions
% cos(kappa_n(z+h)/cos(kappa_n h) where kappa_n are the roots of the 
% dispersion equation for a floating plate
%
%
% the boundary value problem is described in 
% http://www.wikiwaves.org/index.php/Eigenfunction_Matching_for_a_Finite_Floating_Elastic_Plate_using_Symmetry


if nargin == 7
    nu = 0.3;
end

% first we calculate the roots of the dispersion equation for a free
% surface. N is the number of real roots.
k = dispersion_free_surface(alpha,N,h); 

% we need to determine the wavenumber in the y direction
k_y = sin(theta)*k(1);

% Then we calculate the roots of the dispersion equation for a floating
% plate
kappa = dispersion_elastic_surface(alpha,beta,gamma,N+2,h);

% we now find the modified coefficients
k_hat = sqrt(k.^2 - k_y^2);
kappa_hat = sqrt(kappa.^2 - k_y^2);

% we calculate the matrix of inner products. 
A = zeros(N+1,N+1);
B = zeros(N+1,N+3);

for j = 1:N+1;
    A(j,j) = innerproduct_h_to_0(k(j),k(j),h)/cos(k(j)*h)^2;
	for l = 1:N+3; 
		B(j,l) = innerproduct_h_to_0(k(j),kappa(l),h)/...
                 (cos(kappa(l)*h)*cos(k(j)*h));
	end
end


% We define the different blocks of the matrix sym first
M1 = A;
M2 = -B;
M3 = zeros(2,N+1);
M4 = zeros(2,N+3);
for j = 1:N+3
    M4(1,j) = kappa(j)*(kappa_hat(j)^3-kappa_hat(j)*k_y^2*(2-nu))*tan(kappa(j)*h)*tanh(kappa_hat(j)*L);
    M4(2,j) = kappa(j)*(kappa_hat(j)^2-k_y^2*nu)*tan(kappa(j)*h);
end
M5 = zeros(N+1,N+1);
M6 = zeros(N+1,N+3);
for j = 1:N+1;
	for l = 1:N+3; 
		M6(j,l) = (k_hat(j) + tanh(kappa(l)*L)*kappa_hat(l))*B(j,l);
	end
end

% We now assemble the full matrix defining the linear system to solve
M = [M1,M2;M3,M4;M5,M6];

% We can define the right hand term of the equation as well
RHT = zeros(2*N+4,1);
RHT(1) = -A(1,1);
RHT(N+4) = 2*k_hat(1)*A(1,1);

% Finally we can solve the linear system M*S=RHT
S = M\RHT;

a_s = S(1:N+1);
b_s = S(N+2:2*N+4);

% We define the different blocks of the matrix - antisym second. 
M1 = A;
M2 = -B;
M3 = zeros(2,N+1);
M4 = zeros(2,N+3);
for j = 1:N+3
    M4(1,j) = kappa(j)*(kappa_hat(j)^3-kappa_hat(j)*k_y^2*(2-nu))*tan(kappa(j)*h)*coth(kappa_hat(j)*L);
    M4(2,j) = kappa(j)*(kappa_hat(j)^2-k_y^2*nu)*tan(kappa(j)*h);
end
M5 = zeros(N+1,N+1);
M6 = zeros(N+1,N+3);
for j = 1:N+1;
	for l = 1:N+3; 
		M6(j,l) = (k_hat(j) + coth(kappa(l)*L)*kappa_hat(l))*B(j,l);
	end
end

% We now assemble the full matrix defining the linear system to solve
M = [M1,M2;M3,M4;M5,M6];

% We can define the right hand term of the equation as well
RHT = zeros(2*N+4,1);
RHT(1) = -A(1,1);
RHT(N+4) = 2*k_hat(1)*A(1,1);

% Finally we can solve the linear system M*S=RHT
S = M\RHT;

a_a = S(1:N+1);
b_a = S(N+2:2*N+4);

%--------------------------------------------------------------------------

% a set of commands used to test this program.  We plot the potential and
% its derivative on each side and show that they match.  
% uncomment these lines to run the test (except this one)
% z = -h:h/200:0;
% 
% phi_minus = cos(k(1)*(z+h))/cos(k(1)*h); 
% phi_plus = zeros(size(phi_minus));
% 
% for count = 1:N+1
%     phi_minus = phi_minus + a_s(count)*cos(k(count)*(z+h))/cos(k(count)*h);
% end
% for count = 1:N+3
%     phi_plus = phi_plus + b_s(count)*cos(kappa(count)*(z+h))/...
%         cos(kappa(count)*h);
% end
% 
% plot(z,abs(phi_minus),z,abs(phi_plus))
% 
% pause
% 
% phi_minus_n = -k_hat(1)*cos(k(1)*(z+h))/cos(k(1)*h); 
% phi_plus_n = zeros(size(phi_minus));
% 
% for count = 1:N+1
%     phi_minus_n = phi_minus_n + a_s(count)*k_hat(count)*...
%         cos(k(count)*(z+h))/cos(k(count)*h);
% end
% for count = 1:N+3
%     phi_plus_n = phi_plus_n - b_s(count)*kappa_hat(count)*tanh(kappa_hat(count)*L)*...
%         cos(kappa(count)*(z+h))/cos(kappa(count)*h);
% end
% 
% plot(z,abs(phi_minus_n),z,abs(phi_plus_n))
% 
% pause
% 
% phi_minus = cos(k(1)*(z+h))/cos(k(1)*h); 
% phi_plus = zeros(size(phi_minus));
% 
% for count = 1:N+1
%     phi_minus = phi_minus + a_a(count)*cos(k(count)*(z+h))/cos(k(count)*h);
% end
% for count = 1:N+3
%     phi_plus = phi_plus + b_a(count)*cos(kappa(count)*(z+h))/...
%         cos(kappa(count)*h);
% end
% 
% plot(z,abs(phi_minus),z,abs(phi_plus))
% 
% pause
% 
% phi_minus_n = -k_hat(1)*cos(k(1)*(z+h))/cos(k(1)*h); 
% phi_plus_n = zeros(size(phi_minus));
% 
% for count = 1:N+1
%     phi_minus_n = phi_minus_n + a_a(count)*k_hat(count)*...
%         cos(k(count)*(z+h))/cos(k(count)*h);
% end
% for count = 1:N+3
%     phi_plus_n = phi_plus_n - b_a(count)*kappa_hat(count)*coth(kappa_hat(count)*L)*...
%         cos(kappa(count)*(z+h))/cos(kappa(count)*h);
% end
% 
% plot(z,abs(phi_minus_n),z,abs(phi_plus_n))
% 
% % finally we can check that |R|^2 + |T|^2 = 1
% R = 1/2*(a_s(1) + a_a(1));
% T = 1/2*(a_s(1) -  a_a(1));
% 
% abs(R)^2 + abs(T)^2
%--------------------------------------------------------------------------

% this is a program to calculate the inner product of the eigenfunctions
% with each other. 