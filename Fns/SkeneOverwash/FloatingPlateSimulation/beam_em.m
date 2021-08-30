function out = beam_em(lmb,l,xp)
%-------------------------------------------------------------------------
% Function name: beam_em 
% Description: calculates 2D beam/plate dry natural modes of vibration 
% Note: natural modes are normalized in order to have inner product equal 
%       to Kronecker symbol
%-------------------------------------------------------------------------

% Variables
% lmb   :   mode eigenvalue
% l     :   beam/plate half-length
% xp    :   x coord. (calculation points) vector
% out   =   value of the normalized eigenfunction.


out = zeros(length(lmb),length(xp));

for j=1:length(lmb)
    for k=1:length(xp)
        if j-1 == 0               % 1st (rigid) mode
            out(j,k) = 1/sqrt(2*l);
        elseif j-1 == 1           % 2nd (rigid) mode
            out(j,k) = xp(k)/l*sqrt(3/(2*l));
        elseif rem(j-1,2) == 0    % symmetric modes
            out(j,k) = 1/sqrt(2*l)*(cos(lmb(j)*xp(k))/cos(lmb(j)*l)+cosh(lmb(j)*xp(k))/cosh(lmb(j)*l));
        else                      % skew-symmetric modes
            out(j,k) =1/sqrt(2*l)*(sin(lmb(j)*xp(k))/sin(lmb(j)*l)+sinh(lmb(j)*xp(k))/sinh(lmb(j)*l));
        end
    end
end




