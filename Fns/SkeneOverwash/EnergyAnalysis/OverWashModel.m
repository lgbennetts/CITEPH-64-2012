clear all;
close all;
MatFile_NM = strcat('./Data/OW/TransOWMat4.mat');
load(MatFile_NM);


%Plot steepness in vs steepness out
% figure();
% Tp95ka = Ya; %Yka(:,14);
% TM = TMat;%TMat(:,14);
% TowM = TowMat; %TowMat(:,14);
% plot(Tp95ka,TM .*Tp95ka,'--r' )
% hold on;
% plot(Tp95ka,TowM .*Tp95ka,'-k' )
% xlabel('amplitude in (ka)')
% ylabel('amplitude out (ka)')
% 
% T2Mat = IterateMat(Xtp,Ya,TMat,2);
% Tow2Mat = IterateMat(Xtp,Ya,TowMat,2);
% 
% figure();
% Tp95ka = Ya; %Yka(:,14);
% TM = T2Mat;%TMat(:,14);
% TowM = Tow2Mat; %TowMat(:,14);
% plot(Tp95ka,TM .*Tp95ka,'--r' )
% hold on;
% plot(Tp95ka,TowM .*Tp95ka,'-k' )
% xlabel('amplitude in (ka)')
% ylabel('amplitude out (ka)')


T2Mat = TMat.^2;
Tow2Mat = IterateMat(Xtp,Ya,TowMat,2);
%plot the surface - (steepness , Tp) -> Tow / TMat
figure();
surf(Xtp,Ya,T2Mat)
hold on;
% surf(Xtp,Yka, (T2Mat - Tow2Mat).*Yka)
surf(Xtp,Ya, Tow2Mat)
xlabel('Period (s)')
ylabel('Amplitude (a)')
zlabel('Transmission Coefficient')
cb=colorbar; 
set(get(cb,'ylabel'),'String','Transmission Coefficient'); 
title('Comparison of Overwash and Linear Theory Iter - 2')


function IterMat = IterateMat(X,Y,Mat,Iter)

% IterMat = Mat.*interp2(X,Y,Mat,X,Mat.*Y,'makima',0)
% IterMat = IterMat.*interp2(X,Y,Mat,X,IterMat.*Y,'makima',0)
% IterMat = IterMat.*interp2(X,Y,Mat,X,IterMat.*Y,'makima',0)
% IterMat = IterMat.*interp2(X,Y,Mat,X,IterMat.*Y,'makima',0)

if size(X,2) == 1
    IterMat = Mat;%.*interp2(X,Y,Mat,X,Mat.*Y,'makima',0);
    for i = 2:Iter
        IterMat = IterMat.*interp1(Y,Mat,IterMat.*Y,'makima','extrap');
    end
else

    IterMat = Mat;%.*interp2(X,Y,Mat,X,Mat.*Y,'makima',0);
    for i = 2:Iter
        IterMat = IterMat.*interp2(X,Y,Mat,X,IterMat.*Y,'makima',0);
    end

end

return
end

