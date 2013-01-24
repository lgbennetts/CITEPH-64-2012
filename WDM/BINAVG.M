function  x = binavg( x, n )
%
%  x = binavg( x, n )
%
%  Averages the columns of x over n samples and decimates by a factor of n.
%  If x is a vector, its sense (e.g. row or column) is preserved.

%  E.Terray  (based on a program by W. Drennan)  version: 2 June, 1990

[nr,nc] = size( x ) ;

                             % first consider the case when x is a vector
if( nr==1 | nc==1 )
    nx = fix( max([nr,nc]) / n ) ;     % dimension of x after decimation
    y = zeros(n,nx) ;                  % columns of y hold points to average
    y(:) = x(1:n*nx) ;                 % truncate leftover points in x
    x = mean( y ) ;                    % perform the average
    if( nc == 1 )  x = x(:);  end      % make x a column if initially so

else                         % now consider the case where x is an array
    nx = fix( nr / n ) ;
    y = zeros(n,nx*nc) ;               % columns of y hold points to avg
    x = x(1:n*nx,:) ;                  % truncate leftover points in x
    y(:) = x(:) ;                      % resize x into y
    x = zeros(nx,nc) ;                 % template for x after averaging
    x(:) = mean( y ).' ;               % perform the average
end
