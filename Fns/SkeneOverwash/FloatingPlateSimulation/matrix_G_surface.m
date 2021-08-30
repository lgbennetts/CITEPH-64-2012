function G = matrix_G_surface(alpha,H,L,N)
% G = matrix_G_surface(alpha,N,L,H) 
% where alpha is the frequency squared, H is the water depth, L is half the
% plate lenth and 2N+1 is the number of points. The points are evenly space
% starting with -L and ending with L (note that because we use Simpson's
% rule we require and odd number of points.



%first we do the non-zero values
h=2*L/(2*N);
x= 0:h:2*L;

Greenvalues = [0,two_d_finite_depth_Green_surface(alpha,x(2:2*N+1),H)];
weight = 4*h/3*ones(1,2*N+1);
weight(1) = h/3;weight(2*N+1) = h/3;weight(3:2:2*N-1)= 1/2*weight(3:2:2*N-1);


G = zeros(2*N+1,2*N+1);
logvalues = zeros(2*N+1,2*N+1);

for j=1:2*N+1
    for k = 1:j-1
        G(j,k) = weight(k)*Greenvalues(abs(j-k)+1);
        logvalues(j,k) = weight(k)*log(abs(x(j) - x(k)));
    end
    for k = j+1:2*N+1
        G(j,k) = weight(k)*Greenvalues(abs(j-k)+1);
        logvalues(j,k) = weight(k)*log(abs(x(j) - x(k)));
    end
end

% now we need to calculate the solution for the case when  j = k.
sumlogvalues = sum(logvalues.');

% now we find the value of lim R -> 0 G(alpha,R,H) - 1/pi*log(R)

R1 = 1e-1;R2 = 1e-2;
Gsmall = two_d_finite_depth_Green_surface(alpha,[R1,R2],H);
while abs(Gsmall(1) - 1/pi*log(R1) - Gsmall(2) + 1/pi*log(R2)) > 1e-3
%    Gsmall(1) - 1/pi*log(R1);
%    Gsmall(2) - 1/pi*log(R2);
   R1 = R1/10;
   R2 = R2/10;
   Gsmall = two_d_finite_depth_Green_surface(alpha,[R1,R2],H);
   
end

limitvalue = Gsmall(2) - 1/pi*log(R2);




%now we do the j = k values 
G(1,1) = weight(1) * (limitvalue) + (2/pi)*(L*log(2*L)-L) - 1/pi*sumlogvalues(1);
for j=2:2*N
    G(j,j) = weight(j) * (limitvalue) + (1/pi)*((2*L-x(j))*log(2*L-x(j))- 2*L + (x(j))*log(x(j)))  ...
        - 1/pi*sumlogvalues(j);
end
G(2*N+1,2*N+1) = weight(2*N+1) * (limitvalue) + (2/pi)*(L*log(2*L)-L) - 1/pi*sumlogvalues(2*N+1);