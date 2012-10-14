function L = fn_LogAnalyticp(nn,mm,vars)

% - L = fn_LogAnalytic_pm(nn,mm,vars)
%
% - inner prod wrt exp(i*nn*th-i*mm*tau)
% - vars = [~,width,~,Rad]
% - Tols = [~,~,Rtol]

width = vars(2); % - y \in (-width,width)
Rad = vars(4); % - the radius of the int contour

pp = pi/width; Rp = Rad*pp;

if and(nn==0,mm==0)
 L = pi*log(Rp);   
elseif nn==mm
 L = -pi/2/abs(nn);
else
 L = 0;
end
