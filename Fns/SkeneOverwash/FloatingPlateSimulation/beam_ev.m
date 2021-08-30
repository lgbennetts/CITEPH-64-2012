function lmb = beam_ev(n,l)
%-------------------------------------------------------------------------
% Function name: beam_ev 
% Description: finds eigenvalues for the 2D beam/plate dry vibrations  
% Note: eigenvalues for the rigid modes are also included
%-------------------------------------------------------------------------

% Variables
% n     :   number of the desired eigenvalues
% l     :   beam/plate half-length 
% lmb   =   vector of eigenvalues

lmb = zeros(1,n);
p=0;
q=0;
for k=1:n
    if rem(k-1,2) == 0
%         temp = fzero(@(x) tan(x)+tanh(x),p*pi)
%         lmb(k) = temp/l;
%         p=p+1; 
            funct=@(x) tan(x)+tanh(x);
            lmb(k)=goldensectionsearch(funct,(p-1)*pi+pi/2,(p)*pi+pi/2)/l;
            p=p+1;
    else
%         temp = fzero(@(x) -1*tan(x)+tanh(x),q*pi)
%         lmb(k) = temp/l;
%         q=q+1;
            funct=@(x) -1*tan(x)+tanh(x);
            lmb(k)=goldensectionsearch(funct,(q-1)*pi+pi/2,(q)*pi+pi/2)/l;
            q=q+1;
    end
end
end

function [value] = goldensectionsearch(func,xL,xR)
    golden=2/(1+sqrt(5));
    x2=xL+(xR-xL)*golden;
    x1=xR-(xR-xL)*golden;
    func1=func(x1);
    func2=func(x2);
    threshold=0.0001;
    while (xR-xL)>threshold
        if (func1>func2 && func2>0) || (func2>func1 && func2<0)
            xL=x1;
            x1=x2;
            x2=xL+(xR-xL)*golden;
            func1=func2;
            func2=func(x2);
        else
            xR=x2;
            x2=x1;
            x1=xR-(xR-xL)*golden;
            func2=func1;
            func1=func(x1);
        end
    end
    value=(xR+xL)/2;
end