function theta = GrafThetaVal(v, u, delta)

% - function to calculate the value of theta_pm,pm:
%                                   
%                               tan(theta_pm,+)= pm v/u
%                               cos(theta_pm,+)= u/k
%                               sin(theta_pm,+)= pm v/k
%   and
%                                   
%                               tan(theta_pm,-)= mp v/u
%                               cos(theta_pm,-)= -u/k
%                               sin(theta_pm,-)= pm v/k

% - nb. k need not be known here - %

% - del is a vector length 2 of entries +1 or -1 - %

% - fn arctan will retrieve root in (-pi/2,pi/2) - %

if u == 0 % case u=0 -> v > 0 in R
    
    theta = delta(1)*pi/2;
    
else

% - case 1: v \in Real

if imag(v) == 0

if delta(2)*u > 0 % - u can be either +ve or -ve  
    
    theta = atan(delta(2)*delta(1)*v/u);
    
elseif delta(1) > 0 % - nb v must be +ve
    
    theta = atan(delta(2)*v/u);
    
    theta = theta + pi;
    
elseif delta(1) < 0
    
    theta = atan(-1*delta(2)*v/u);
    
    theta = theta - pi;
    
end

% - case 2 v \in i*Real

else
    
    % - no periodicity here
    
    if delta(2)*u > 0
    
    theta = atan(delta(2)*delta(1)*v/u);
    
    elseif delta(2)*u < 0
        
    theta = pi + atan(delta(2)*delta(1)*v/u);
    
    end
    
end

end
    
    
    
    
    