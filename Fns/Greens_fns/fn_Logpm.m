function Logpm = fn_Logpm(xm,ypm,vars,Tols)

% - xm = x-x0 
%   yp = y \pm y0
% - vars = [~,width,~]

width = vars(2); % - y \in (-width,width)
% kk = vars(3); % - wavenumber
pp = pi/width;

Ep = exp(-pp*abs(xm)+1i*pp*ypm);
Em = exp(-pp*abs(xm)-1i*pp*ypm);

if (1-Em)*(1-Ep)<Tols(3)
    Logpm=0;
else
    Logpm = log((1-Em)*(1-Ep))/4/pi;
end

%Logpm = 1i*log((1-Em)*(1-Ep))/pp;
%Logpm = Logpm/8i/width;

return