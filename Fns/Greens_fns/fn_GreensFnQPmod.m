function [G, g_sk] = fn_GreensFnQPmod(xx,x0,vars,skip,Tols)

% - 27.05.10 - after noting mistake in Green's fn

e0 = Term0(xx(1)-x0(1),vars);

[Gm, g_sk] = fn_Gpm_til(xx(1)-x0(1),xx(2)-x0(2),vars,skip,Tols);

Logm = fn_Logpm(xx(1)-x0(1),xx(2)-x0(2),vars,Tols);

G = (e0 + Gm) + Logm;  

return

%% - Subfns - %%

function ee = Term0(xm,vars)

width = vars(2); kk = vars(3);

ee = exp(1i*kk*abs(xm))/4i/width/kk;

return