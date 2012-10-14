function Logpm = fn_LogSing(xx,x0,c0,pm,vars,Tols)

% 21.09.12 FOR TESTING ONLY
%
% - xx = (x,y) -> field vars
%   x0 = (x0,y0) -> source vars
% - c0 = (x,y) is the centre of the circle
% - pm = +1 y+y0, -1 y-y0
% - vars = [~,width,~,Rad]
% - Tols = [~,~,Rtol]

width = vars(2); % - y \in (-width,width)
pp = pi/width;
Rad = vars(4); % - the radius of the int contour

xm = xx(1)-x0(1); ypm = xx(2)+pm*x0(2);

th = angle(xx(1)-c0(1)+1i*(xx(2)-c0(2))); 
tau = angle(x0(1)-c0(1)+1i*(x0(2)-c0(2)));

gg = (2*(Rad*pp)^2)*(1-cos(th+pm*tau));

if abs(xm+1i*ypm)<Tols(3)
    Logpm = 0;
else
    Logpm = log(gg)/4/pi;
end

return