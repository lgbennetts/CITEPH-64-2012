function G = fn_GreensFnQuasiPer_Basic(xx,x0,vars,Tols)

% - FOR TESTING ONLY

e0 = Term0(xx(1)-x0(1),vars);

Gm = fn_Gp_Chan(xx(1)-x0(1),xx(2)-x0(2),vars,Tols);

G = e0 + Gm; 

return

%% - Subfns - %%

function ee = Term0(xm,vars)

width = vars(2); kk = vars(3);

ee = exp(1i*kk*abs(xm))/4i/width/kk;

return