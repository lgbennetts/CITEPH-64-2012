function out = two_d_finite_depth_Green_surface(alpha,R,H)
% out = two_d_finite_depth_Green_surface(alpha,R,H)
% calculates the finitedepth green function at the for both source and 
% field point on the free surface for two dimensions.   
% R is the distance between two points and H is the water depth. 
% The values of R must
% be ordered from the smallest to the largest and be positive.
%
% Details can be found on 
% http://www.wikiwaves.org/index.php/Free-Surface_Green_Function

error = 1e-10; % this is the difference between the actual and
% estimated answer which we are prepared to tolerate. 

%first of all we calculate the roots of the dispersion equation 

N=10;

out = zeros(size(R));

mroots = dispersion_free_surface(alpha,N,H);



for j = 1:length(R)
         out(j) =  -sum((exp(-mroots*R(j)))./(tan(mroots*H) + H*mroots.*sec(mroots*H).^2));
        % Now we determine if the series has converged. 
        last_term =((exp(-mroots(N+1)*R(j)))./(tan(mroots(N+1)*H) + H*mroots(N+1).*sec(mroots(N+1)*H).^2));
         %abs(last_term / out(j))
        while abs(last_term / out(j)) > error
           %we have not converged
           N = N*2;
           mroots = dispersion_free_surface(alpha,N,H);
            out(j) =  -sum((exp(-mroots*R(j)))./(tan(mroots*H) + H*mroots.*sec(mroots*H).^2));
        % Now we determine if the series has converged. 
        last_term =((exp(-mroots(N+1)*R(j)))./(tan(mroots(N+1)*H) + H*mroots(N+1).*sec(mroots(N+1)*H).^2));
        end

           
     
end