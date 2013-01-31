function [Gpm,gp] = fn_Gpm_til(xm,ypm,vars,skip,Tols)

% - With the logarithms removed
% -> These are replaced by the functions log_pm
%
% - xm = x-x0 
%   ypm = y \pm y0
% - vars = [~,width,kk]
% - Tols = [MaxAz,Convergence Stop]

% - 27.05.10 - after noting mistake in Green's fn
% - 26.04.10 - adjusted for resonances -> nb. always 2 as normal incidence

% - last modified: 31.01.13

%USE_POLYLOGS = 1;

width = vars(2); % - y \in (-width,width)
kk = vars(3); % - wavenumber
pp = pi/width;

Gpm = 0; 

loop=1; %enough=0; sv_vec = 8*width*Tols(2)*ones(3,1);

while loop<=Tols(1) %and(enough == 0, loop<=Tols(1))
    
    mu_m = loop*pp;
    u_m = sqrt(kk^2-mu_m^2);
    
    sm = exp(1i*u_m*abs(xm))*(...
          exp(1i*mu_m*ypm) + exp(-1i*mu_m*ypm) );
      
    %%% REMOVE RESONANT TERM %%% 
    if loop==skip
     sm=0;
    else
     sm = sm/u_m; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tm = exp(-mu_m*abs(xm))*(...
          exp(1i*mu_m*ypm) + exp(-1i*mu_m*ypm) );
      
    tm = tm/mu_m/1i;
    
    Gpm = Gpm + sm-tm;
    
%     sv_vec(1:2) = sv_vec(2:3); sv_vec(3) = abs(sm-tm); 
%     
%     if max(sv_vec)<Tols(2)
%       enough = 1;
%     end
            
    loop = loop + 1;
    
end

% if loop>Tols(1)
%     display(['terms=',num2str(loop)])
% end

Gpm = Gpm/4i/width; 

%% EXTRA INFO IF RESONANCE OCCURS

if isempty(skip)==1
 gp=[];
else
 %%% THE MISSING TERM *u_m
 mu_m = skip*pp;
 u_m = sqrt(kk^2-mu_m^2);
    
 sm = exp(1i*u_m*abs(xm))*(...
  exp(1i*mu_m*ypm) + exp(-1i*mu_m*ypm) );
      
 gp = sm/4i/width;
end

return    