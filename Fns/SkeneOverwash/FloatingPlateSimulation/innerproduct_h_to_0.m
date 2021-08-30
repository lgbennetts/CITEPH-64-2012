function out = innerproduct_h_to_0(k,kappa,h)
% program to calculate 
% \int_{-h}^{0}\cos \kappa\left( z+h\right) 
% \cos k\left(z+h\right) dz

% first we test for the case when \kappa and k are equal (or nearly equal)
if abs((k-kappa)/h) > 1e-6
%     out = 1/2 * (sin(- kappa*h + k*h)) / (k - kappa) ...
%         + 1/2 * (sin(kappa*h + k*h)) / (kappa + k);
      out = (k*sin(k*h)*cos(kappa*h)-kappa*cos(k*h)*sin(kappa*h))/...
          (k^2-kappa^2);
else
    if k == 0
        out = h-d;
    else
        out = (1/2)*((cos(k*h)*sin(k*h) + k*h)/k);
    end
end

