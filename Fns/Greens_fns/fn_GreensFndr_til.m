function Gr = fn_GreensFndr_til(xx,x0,c0,vars,Tols)

Gm_dr= fn_Gpm_dr_til(xx,x0,c0,vars,Tols);

Logm_dr = fn_Logpm_dr(xx,x0,c0,vars,Tols);

Gr = Gm_dr + Logm_dr;

return