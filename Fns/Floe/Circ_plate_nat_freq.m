% FABIEN'S CODE

function [freq,shape] = Circ_plate_nat_freq(nu,N,J)

if nargin == 0
    N = 10;
    J = 5;
    nu = 0.33;
elseif nargin == 1
    N = 10;
    J = 5;
end

method = 'bisection';

freq = zeros(N+1,J+1);
for n = 0:N
    freq(n+1,:) = FindRoots(n,nu,J,method);
end

shape = zeros(N+1,J+1);
for n = 1:N+1
    for j = 1:J+1
        shape(n,j) = F1(freq(n,j),n-1,nu)/F2(freq(n,j),n-1,nu);
    end
end

return

function mroot = FindRoots(n,nu,J,method)
%calculates the jth root of the nth angular mode using the homotopy method

mroot = zeros(1,J+1);
for j =0:J
    if j == 0;
        if n>1
            if strcmp(method,'newton')
                mroot(j+1) = newton(n,nu,n);
            elseif strcmp(method,'secant')
                mroot(j+1) = secant(n,nu,n,pi*n/2);
            elseif strcmp(method,'bisection')
                mroot(j+1) = bisection(n,nu,n,pi*n/2);
            end
        end
    else
        if strcmp(method,'newton')
            mroot(j+1) = newton(n,nu,pi/2*(n+2*j));
        elseif strcmp(method,'secant')
            mroot(j+1) = secant(n,nu,pi/2*(n+2*(j-1)),pi/2*(n+2*j));    
        elseif strcmp(method,'bisection')
            mroot(j+1) = bisection(n,nu,pi/2*(n+2*(j-1)),pi/2*(n+2*j));
        end
    end
end

function out = bisection(n,nu,guess_min,guess_max)

if ~exist('COMM','var'); COMM=0; end

if guess_min == 0
    guess_min = guess_min+0.01;
end
if COMM
 if f(guess_min,n,nu)*f(guess_max,n,nu)>0
    fprintf('no roots found on this interval')
 end
end
ans1 = guess_min;
ans2 = guess_max;
while abs(ans2 - ans1) > 1e-9
    fval = f((ans1+ans2)/2,n,nu);
    if fval*f(ans1,n,nu) < 0
        ans2 = (ans1+ans2)/2;
    else
        ans1 = (ans1+ans2)/2;
    end
end
out = ans2;

function out = secant(n,nu,guess_min,guess_max)
if guess_min == 0
    guess_min = guess_min+0.01;
end
ans2 = guess_min;
out = guess_max;
while abs(ans2 - out) > 1e-9
    ans1 = ans2;
    ans2 = out;
    out = ans2 - (ans2 - ans1)*f(ans2,n,nu)/(f(ans2,n,nu)-f(ans1,n,nu));
end

function out = newton(n,nu,guess)
%calculates the root nearest the root guess.
ans1 = guess+1;
out = guess;
while abs(ans1 - out) > 1e-9
    ans1 = out;
    out = ans1 - f(ans1,n,nu)/df(ans1,n,nu);
%     if out > pi/2*(n+2*j+0.5)
%         out = ans1-0.5;
%     end
end


function out = f(lambda,n,nu)
% out = F1(lambda,n,nu)/F2(lambda,n,nu) - G1(lambda,n,nu)/G2(lambda,n,nu);
out = F1(lambda,n,nu)*G2(lambda,n,nu) - G1(lambda,n,nu)*F2(lambda,n,nu);

function out = df(lambda,n,nu)
out = (DF1(lambda,n,nu)*F2(lambda,n,nu) - ...
    F1(lambda,n,nu)*DF2(lambda,n,nu))/F2(lambda,n,nu)^2 - ...
    (DG1(lambda,n,nu)*G2(lambda,n,nu) - ...
    G1(lambda,n,nu)*DG2(lambda,n,nu))/G2(lambda,n,nu)^2;

function out = F1(lambda,n,nu)
out = lambda.^2.*besselj(n,lambda)+(1-nu)*(lambda.*dbesselj(n,lambda) - ...
    n^2*besselj(n,lambda));

function out = F2(lambda,n,nu)
out = lambda.^2.*besseli(n,lambda)-(1-nu)*(lambda.*dbesseli(n,lambda) - ...
    n^2*besseli(n,lambda));

function out = G1(lambda,n,nu)
out = lambda.^3.*dbesselj(n,lambda)+...
    n^2*(1-nu)*(lambda.*dbesselj(n,lambda)-besselj(n,lambda));

function out = G2(lambda,n,nu)
out = lambda.^3.*dbesseli(n,lambda)-...
    n^2*(1-nu)*(lambda.*dbesseli(n,lambda)-besseli(n,lambda));

function out = DF1(lambda,n,nu)
out = lambda^2*dbesselj(n,lambda) + 2*lambda*besselj(n,lambda) + ...
    (1-nu)*(lambda*ddbesselj(n,lambda) + dbesselj(n,lambda) - ...
    n^2*dbesselj(n,lambda));

function out = DF2(lambda,n,nu)
out = lambda^2*dbesseli(n,lambda) + 2*lambda*besseli(n,lambda) - ...
    (1-nu)*(lambda*ddbesseli(n,lambda) + dbesseli(n,lambda) - ...
    n^2*dbesseli(n,lambda));

function out = DG1(lambda,n,nu)
out = lambda^3*ddbesselj(n,lambda) + 3*lambda^2*dbesselj(n,lambda) + ...
    n^2*(1-nu)*lambda*ddbesselj(n,lambda);

function out = DG2(lambda,n,nu)
out = lambda^3*ddbesseli(n,lambda) + 3*lambda^2*dbesseli(n,lambda) - ...
    n^2*(1-nu)*lambda*ddbesseli(n,lambda);

% ADDED BY LUKE

function out = dbesselj(n,z,scal)
% this function calcualtes the derivative of besselj(n,z)
if nargin==3
    out = (besselj(n-1,z,scal) - besselj(n+1,z,scal))/2;
else
    out = (besselj(n-1,z) - besselj(n+1,z))/2;
end

function out = dbesseli(n,z,scal)
% this function calcualtes the derivative of besseli(n,z)
if nargin==3
    out = (besseli(n+1,z,scal) + besseli(n-1,z,scal))/2;
else
    out = (besseli(n+1,z) + besseli(n-1,z))/2;
end