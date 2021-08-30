clear all;
close all;



%Models - Transmission Single Floe
th_res=100;
SURGE=1;
terms_grn=100;
extra_pts=[];

if ~exist('Param','var'); Param = ParamDef_Oceanide(4); 
Param = ModParam_def(Param,1,1e2,extra_pts,terms_grn,th_res); end

%Linear Theories

dum_c = pi*((Param.floe_diam/2)^2)/2;
EDN_Real = (dum_c*Param.MIZ_length)/Param.floe_diam;
EDN_2D = (dum_c*Param.MIZ_length)/Param.floe_diam/2;


MatFile_NM = strcat('./Data/OW/TransOWMat.mat');
load(MatFile_NM);

SingYa = Ya;
SingTM = TMat;
SingTMOW = TowMat;

Tp = Xtp(1,1);
freq=1/Tp; %wave frequency Hz
omega=2*pi*freq; %angular frequency rad/s
alpha=omega^2/Param.g; %scaled frequency
k= dispersion_free_surface(alpha,0,Param.bed); %wavenumber
kr = imag(k);

%Plot steepness in vs steepness out
T2Mat = SingTM.^EDN_Real;
Tow2Mat = IterateMat(Xtp,Ya,SingTMOW,2);

TransDataReg = load('./Data/Gen/NewATrans39Reg');
IndexMatchTp = find([TransDataReg.TpA] == Tp) ; 

DataIn = TransDataReg.CalTra(IndexMatchTp);
DataOut = TransDataReg.ExpTra(IndexMatchTp);

figure();
set(gca,'FontSize',24) 
plot(SingYa,T2Mat .*SingYa,'--r','DisplayName', 'Linear Model' )
hold on;
plot(SingYa,Tow2Mat .*SingYa,'-k','DisplayName', 'Overwash Model' )
plot(DataIn,DataOut,'ob','DisplayName', 'Regular Data')
xlabel('amplitude in (m)')
ylabel('amplitude out (m)')
title(['Tp = ', num2str(Tp), ' cL/a: ', num2str(EDN_Real)]);
legend()


% T2Mat = IterateMat(Xtp,Yka,TMat,2);
% Tow2Mat = IterateMat(Xtp,Yka,TowMat,2);
% %plot the surface - (steepness , Tp) -> Tow / TMat
% figure();
% % surf(Xtp,Yka,T2Mat)
% hold on;
% % surf(Xtp,Yka,TMat.^2)
% surf(Xtp,Yka, (T2Mat - Tow2Mat).*Yka)
% % surf(Xtp,Yka, TowMat.^2)
% xlabel('Period (s)')
% ylabel('Steepness (ka)')
% zlabel('Transmission Coefficient')
% cb=colorbar; 
% set(get(cb,'ylabel'),'String','Transmission Coefficient'); 
% title('Comparison of Overwash and Linear Theory Iter - 2')


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

