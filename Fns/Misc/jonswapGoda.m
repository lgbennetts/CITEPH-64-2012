function [y] = jonswapGoda(ff,fp)


alfa=0.0081;
gamma=3.3;
g= 9.81;

ffs = ff./fp;
sigma       = (ffs < 1) * 0.07 + (ffs >= 1) * 0.09;
Term1 = alfa*g^2*(2*pi)^(-4)*ff.^(-5);
Term2 = exp(-5/4*(ffs).^(-4));
Term3 = gamma.^(exp((-1/2)*(ffs - 1)./sigma).^2);


y = Term1.*Term2.*Term3;

% %% create output
% if y(end)/max(y)>.01
% 	warning(sprintf('%s\n%s','Upper limit domain is too small to cover full spectrum.',...);
% 	'Consider raising upper limit.')); 
% end 
% if ~OPTset.wp & ~OPTset.Hs
% 	y = y / max(y);
% elseif OPTset.wp & OPTset.Hs
% 	y = y/trapz(x,y)*OPT.Hs^2*pi/8; 
% else
% 	y = y;
% end

