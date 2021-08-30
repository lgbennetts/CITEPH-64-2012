
close all;

Unit = 0:0.2:1;

[X,Y] = meshgrid(Unit);

Z = 0.5*X;
figure();
surf(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('f(x,y)')

figure();
Z1 =  interp2(X,Y,Z,X,Z.*Y,'makima',0);
% Z2 = interp2(X,Y,Z1,X,Z1.*Y,'makima',0);
surf(X,Y,Z1)
xlabel('x')
ylabel('y')
zlabel('f(x,f(x,y))')
title('f(x,f(x,y))')



% function IterMat = IterateMat(X,Y,Mat,Iter)
% 
% % Xo = X(1,:);
% % Yo = Y(:,1);
% % InterExtrap = [0:0.01:2,OWData.kI(2:end)];
% % InterExtrap = sort(InterExtrap);
% IterMat =  interp2(X,Y,Mat,X,Y,'makima',0);
% 
% for i = 2:Iter
%     IterMat.*Y;
%     IterMat = interp2(X,Y,Mat,X,IterMat.*Y,'makima',0);
% end
% 
% return
% end

