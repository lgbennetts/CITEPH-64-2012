function k = wavek( f, H, tol )
%
%    k = wavek( f, H, tol )
%
%    Returns a column vector of wavenumbers [k] (1/m) at frequencies [f] (Hz)
%    using the gravity wave dispersion relation for water depth H (m).  tol
%    is an optional error bound on the predicted f's (the default is 1e-4). 
%    Note: [f] > 0 must increase monotonically from small to large.  This
%    routine is optimized for speed and works best when length(f) is large.

                             % some initial housekeeping

if nargin<3,  tol = 1e-4;  end         % default tolerance
g = 9.80171;   c = 2*pi*sqrt(H/g);     % g at 40N (m/s^2)
f = c * f(:);  Nf = length(f);         % non-dimensionalize f
k = zeros(Nf,1);
                             % use approximations for large & small f

f_r = 1.1*sqrt( abs(log(tol/2))/2 );   % define f LARGE if error > 110% * tol
nr = find( f >= f_r ) ;                % if so, use the limiting deep-water
if isempty(nr)                         % dispersion relation  k = f^2
    nr = Nf+1 ;
else
    k(nr) = ( f(nr).^2 ) ;
    nr = min( nr ) ;
end

f_l = 0.7*( tol * 604800/281 )^(1/8);  % f is SMALL if error > 70% * tol
nl = find( f <= f_l ) ;                % if so, use shallow water expansion
if isempty(nl)                         % k = f*(1+f^2/6+11f^4/360+17f^6/5040)
    nl = 0 ;
else
    k(nl) = f(nl) .* polyval([17/5040,11/360,1/6,1], f(nl).^2) ;
    nl = max( nl ) ;
end
                             % calculate the nonlinear k(f) regime by solving
n = (nl+1):1:(nr-1) ;        % f^2 - k*tanh(k) == 0 using Newton-Raphson
if ~isempty( n )
    dk = 2*tol ;
    k(n) = f(n).^2 ;                   % initial guess is deep water limit

    while max( abs(dk) ) > tol         % Newton-Raphson iteration loop
        t = tanh( k(n) ) ;
        dk = -(f(n).^2 - k(n).*t) ./ ( t + k(n).*( 1 - t.^2 ) ) ;
        k(n) = k(n) - dk ;
    end
    k(n) = abs( k(n) ) ;               % f(k) = f(-k), so k>0 == k<0 roots
end

N = min( find( k>50 ) ) ;              % inform user if err > tol
if isempty(N),  N = Nf;  end
err = abs( f(2:N) - sqrt( k(2:N).*tanh(k(2:N)) ) )./f(2:N) ;
if max( err ) > tol,  fprintf('\n WAVEK: error exceeds %g \n', tol),  end

k = k / H ;                            % give k dimensions of 1/meter
