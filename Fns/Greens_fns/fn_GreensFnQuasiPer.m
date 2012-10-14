function [G, g_sk] = fn_GreensFnQuasiPer(xx,x0,c0,vars,skip,Tols)

% - 27.05.10 - after noting mistake in Green's fn

e0 = Term0(xx(1)-x0(1),vars);

[Gm,g_sk] = fn_Gpm_til(xx(1)-x0(1),xx(2)-x0(2),vars,skip,Tols);

Logm_v1 = fn_Logpm_til(xx,x0,c0,-1,vars,Tols);

G = e0 + Gm + Logm_v1; % 

return

%% - Subfns - %%

function ee = Term0(xm,vars)

width = vars(2); kk = vars(3);

ee = exp(1i*kk*abs(xm))/4i/width/kk;

return