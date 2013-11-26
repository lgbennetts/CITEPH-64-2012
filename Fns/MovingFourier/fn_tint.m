function [t0,t1] = fn_tint(Tp,T_pers,t0,tvec)

% if Tp==.65
%  t0=t0+T_pers*Tp;
%  t1=min(tvec)-T_pers*Tp;
% elseif Tp==.8
%  t0=t0+T_pers*Tp;
%  t1=min(tvec)-T_pers*Tp;
% elseif Tp==.95
%  p =0.75;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp; 
% elseif Tp==1.1
%  p =0.75;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% elseif Tp==1.25
%  p =0.5;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% elseif Tp==1.4
%  p =0.5;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% elseif Tp==1.55
%  p =0.5;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% elseif Tp==1.7
%  p =0.5;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% elseif Tp==1.85
%  p =0.5;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% elseif Tp==2
%  p =0.5;
%  t0=t0+p*T_pers*Tp;
%  t1=min(tvec)-p*T_pers*Tp;
% else
%  t0=t0; t1=min(tvec);
% end

t1 = min(tvec); clear tvec

tw = (t1-t0)/3;
tm = (t0+t1)/2;

t0 = tm - tw;
t1 = tm + tw;

return