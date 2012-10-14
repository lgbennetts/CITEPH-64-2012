function Gpm = fn_Gp_Chan(xm,ypm,vars,Tols)

% FOR TESTING

% - xm = x-x0 
%   ypm = y \pm y0
% - vars = [~,width,kk]
% - Tols = [MaxAz,Convergence Stop]

% - last modified: 20.09.12

width = vars(2); % - y \in (-width,width)
kk = vars(3); % - wavenumber
pp = pi/width;

Gpm = 0; 

loop=1; enough=0; sv_vec = 8*width*Tols(2)*ones(3,1);

while and(enough == 0, loop<=Tols(1))
    
    mu_m = loop*pp;
    u_m = sqrt(kk^2-mu_m^2);
    
    sm = exp(1i*u_m*abs(xm))*(...
          exp(1i*mu_m*ypm) + exp(-1i*mu_m*ypm) );
      
    sm = sm/u_m; 
    
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

Gpm = Gpm/4i/width; 

return
    
    