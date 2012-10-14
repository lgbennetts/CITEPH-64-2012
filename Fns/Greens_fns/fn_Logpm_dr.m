function Logpm_dr = fn_Logpm_dr(xx,x0,c0,vars,Tols)

% - Logpm_dr = fn_Logpm_dr(xx,x0,c0,pm,vars,Tols)
%
% - xx = (x,y) -> field vars
%   x0 = (x0,y0) -> source vars
% - c0 = (x,y) is the centre of the circle
% - vars = [~,width,~,Rad]
% - Tols = [~,~,Rtol]

width = vars(2); % - y \in (-width,width)
% kk = vars(3); % - wavenumber
pp = pi/width; Rp = vars(4)*pp;

xm = xx(1)-x0(1); ypm = xx(2)-x0(2);

if xm>=0 % - arbitrarily assign this bias
 sgnx = 1;
else
 sgnx=-1;
end

th = angle(xx(1)-c0(1)+1i*(xx(2)-c0(2))); 

Ep = exp(-pp*abs(xm)+1i*pp*ypm);
Em = exp(-pp*abs(xm)-1i*pp*ypm);

Ethp = exp(1i*sgnx*th); Ethm = exp(-1i*sgnx*th); % - nb. sign(xm)=0 iff xm=0

if abs(xm+1i*ypm)<Tols(3)
    Logpm_dr = (1-Rp*sgnx*cos(th))/Rp; % - nb. sign(xm)=0 iff xm=0
else
    Logpm_dr = sgnx*(Ethm*Ep*(1-Em) + Ethp*Em*(1-Ep))/(1-Ep)/(1-Em); % - nb. sign(xm)=0 iff xm=0
end

Logpm_dr = Logpm_dr/4/width;

return