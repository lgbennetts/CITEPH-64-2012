function Besdz = Bessel_dz(BesselFn, N, z, scale, hank, vargin)

% - derivatives of Bessel functs. using recurrence formulae - %

if nargin == 3

if N == 0
    
    Besdz = - feval( BesselFn, 1, z);
    
else
    
    Besdz = feval( BesselFn, N-1, z) - feval( BesselFn, N+1, z);
    
    Besdz = Besdz / 2;
    
end

elseif nargin == 4
    
if N == 0
    
    Besdz = - feval( BesselFn, 1, z, scale);
    
else
    
    Besdz = feval( BesselFn, N-1, z, scale) - feval( BesselFn, N+1, z, scale);
    
    Besdz = Besdz / 2;
    
end 

elseif nargin == 5
    
if N == 0
    
    Besdz = - feval( BesselFn, 1, hank, z, scale);
    
else
    
    Besdz = feval( BesselFn, N-1, hank, z, scale) - feval( BesselFn, N+1, hank, z, scale);
    
    Besdz = Besdz / 2;
    
end     

end