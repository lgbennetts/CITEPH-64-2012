function x = mod360(x)
%
%  function x = mod360(x)
%
%    Given array of angles x [deg], returns x in range 0-360. i.e. mod(x,360).
%
x=rem(x,360);
x(x<0) = x(x<0) + 360;
x(x>=360) = x(x>=360) - 360;
