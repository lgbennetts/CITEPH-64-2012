function [ xi, eta ] = BendShear3dAxi_fn( parameter_vector, Bes, arg )

p = 1-parameter_vector(3);

% - Bes = [ fn, N ]

% - fn: 1 = J; 2 = Y; 3 = H1; 4 = I

if Bes(1) == 1
    
    if length(Bes)==3
     B = besselj(Bes(2), arg, Bes(3));
     Bz = Bessel_dz(@besselj, Bes(2), arg, Bes(3));
    else
     B = besselj(Bes(2), arg);
     Bz = Bessel_dz(@besselj, Bes(2), arg);
    end

 xi = ( arg.^2 - p*(Bes(2).^2) ).*B + p.*arg.*Bz;
 eta = ( arg.^2 + p*(Bes(2).^2) ).*arg.*Bz - (Bes(2).^2).*p.*B;

% elseif Bes(1) == 2
%     
%     B = bessely(Bes(2), arg);
%     
%     Bz = Bessel_dz(@bessely, Bes(2), arg);
%     
% xi = ( arg.^2 - p.*(Bes(2).^2) ).*B + p.*arg.*Bz;
% eta = ( arg.^2 + p.*(Bes(2).^2) ).*arg.*Bz - (Bes(2).^2).*p.*B;

elseif Bes(1) == 3 % - Hankel of 1st kind !
    
    if length(Bes)==3
     B = besselh(Bes(2), 1, arg, Bes(3));
     Bz = Bessel_dz(@besselh, Bes(2), arg, 1, Bes(3));
    else
     B = besselh(Bes(2), arg);
     Bz = Bessel_dz(@besselh, Bes(2), arg);
    end
    
xi = ( arg.^2 - p.*(Bes(2).^2) ).*B + p.*arg.*Bz;
eta = ( arg.^2 + p.*(Bes(2).^2) ).*arg.*Bz - (Bes(2).^2).*p.*B;

elseif Bes(1) == 4
    
    B = besseli(Bes(2), arg);
    
    Bz = Bessel_Idz(@besseli, Bes(2), arg);
    
xi = ( arg.^2 + p.*(Bes(2).^2) ).*B - p.*arg.*Bz;
eta = ( arg.^2 - p.*(Bes(2).^2) ).*arg.*Bz + (Bes(2).^2).*p.*B;

end