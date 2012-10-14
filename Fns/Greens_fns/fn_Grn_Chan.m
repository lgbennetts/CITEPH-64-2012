function Gpm = fn_Grn_Chan(xm,y,y0,vars,Tols)

% - TESTING ONLY (FIND NOTATION IN DOC Wavetank.pdf (2012))

% - xm = x-x0 
% - vars = [~,width,kk]
% - Tols = [MaxAz,Convergence Stop]

% - last modified: 20.09.12

width = vars(2); % - y \in (-width,width)
kk = vars(3); % - wavenumber
pp = pi/width;

gam_m = 0;
al_m = sqrt(kk^2-gam_m^2);
    
sm = exp(1i*al_m*abs(xm))*cos(gam_m*y)*cos(gam_m*y0);

Gpm = sm/al_m; 

loop=1; enough=0; sv_vec = 8*width*Tols(2)*ones(3,1);

while and(enough == 0, loop<=Tols(1))
    
    gam_m = loop*pp;
    al_m = sqrt(kk^2-gam_m^2);
    
    sm = 2*exp(1i*al_m*abs(xm))*cos(gam_m*y)*cos(gam_m*y0);
      
    sm = sm/al_m; 
    
    Gpm = Gpm + sm;
    
    sv_vec(1:2) = sv_vec(2:3); sv_vec(3) = abs(sm); 
    
    if max(sv_vec)<Tols(2)
      enough = 1;
    end
            
    loop = loop + 1;
    
end

if loop>Tols(1)
    display(['terms=',num2str(loop)])
end

Gpm = Gpm/2i/width; 

return