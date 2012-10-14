function Gpm_dr = fn_Gpm_dr(xx,x0,c0,pm,vars,Tols)

% Gpm = fn_Gpm(xm,ypm,vars,Tols)
% 
% - xx = (x,y) -> field vars
%   x0 = (x0,y0) -> source vars
% - c0 = (x,y) is the centre of the circle
% - vars = [~,width,kk]
% - Tols = [MaxAz,Convergence Stop]

width = vars(2); % - y \in (-width,width)
kk = vars(3); % - wavenumber
pp = pi/width;

xm = xx(1)-x0(1); ypm = xx(2)+pm*x0(2);

if xm>0 % - arbitrarily assign this bias
 sgnx = 1;
elseif xm<0
 sgnx=-1;
else
 sgnx=0; % - seems to work better with this!!! (26.03.10) it's like an average!
end

th = angle(xx(1)-c0(1)+1i*(xx(2)-c0(2))); 

Gpm_dr = 0;

loop=1; enough=0; sv_vec = 8*width*Tols(2)*ones(3,1);

while and(enough == 0, loop<=Tols(1))
    
    mu_m = loop*pp;
    u_m = sqrt(kk^2-mu_m^2);
    
    sm_p = (sgnx*cos(th)+(mu_m/u_m)*sin(th))*exp(1i*u_m*abs(xm)+1i*mu_m*ypm);
    sm_m = (sgnx*cos(th)-(mu_m/u_m)*sin(th))*exp(1i*u_m*abs(xm)-1i*mu_m*ypm);
    
    Gpm_dr = Gpm_dr + sm_p + sm_m;
    
    sv_vec(1:2) = sv_vec(2:3); sv_vec(3) = abs(sm_p + sm_m); 
    
    if max(sv_vec)<Tols(2)
      enough = 1;
    end
            
    loop = loop + 1;
    
end

if 0 %loop>Tols(1)
    display(['Out of terms Gpm_dr, ',num2str(pm)])
end

Gpm_dr = Gpm_dr/8/width;

return
    
    