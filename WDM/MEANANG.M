function ma=meanang(A)
%
%  ma = meanang(A)
%
%  Returns the columnwise angular mean of an array of
%  angles, in radians

ma = atan2(mean(sin(A)),mean(cos(A)));
