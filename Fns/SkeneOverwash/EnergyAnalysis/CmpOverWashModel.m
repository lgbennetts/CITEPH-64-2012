clear all;
close all;
MatFile_NM = strcat('./Data/OW/TransOWMatDavid.mat');
DavidDat = load(MatFile_NM);

MatFile_NM = strcat('./Data/OW/TransOWMatLuke.mat');
LukeDat = load(MatFile_NM);


Tp = DavidDat.Xtp(1,1);
freq=1/Tp; %wave frequency Hz
omega=2*pi*freq; %angular frequency rad/s
alpha=omega^2/9.81; %scaled frequency
k= dispersion_free_surface(alpha,0,0.9); %wavenumber
kr = imag(k);
Param.floe_diam = 0.99;
Param.MIZ_length = 5;
dum_c = pi*((Param.floe_diam/2)^2)/2;
EDN_Real = (dum_c*Param.MIZ_length)/Param.floe_diam;
EDN_2D = (dum_c*Param.MIZ_length)/Param.floe_diam/2;

figure();
set(gca,'FontSize',30) 
plot(DavidDat.Ya*kr,DavidDat.TMat.*DavidDat.Ya*kr,'--r','DisplayName','WikiWaves - Tranmission','LineWidth',2)
hold on;
plot(DavidDat.Ya*kr,DavidDat.TowMat.*DavidDat.Ya*kr,'--k','DisplayName','WikiWaves - Tranmission + Overwash' ,'LineWidth',2)

plot(LukeDat.Ya*kr,LukeDat.TMat.*LukeDat.Ya*kr,'-b','DisplayName','Lukes Code - Tranmission','LineWidth',2 )
plot(LukeDat.Ya*kr,LukeDat.TowMat.*LukeDat.Ya*kr,'-g','DisplayName','Lukes Code - Tranmission + Overwash','LineWidth',2 )
xlabel('steepness in (ka)')
ylabel('steepness out (ka)')
legend()
title(['Comparison Of Linear Models Tp = 0.9'])




%Plot steepness in vs steepness out
LBT2Mat = LukeDat.TMat.^EDN_Real;
LBTow2Mat = IterateMat( LukeDat.Xtp, LukeDat.Ya,LukeDat.TowMat,2);
DST2Mat = DavidDat.TMat.^EDN_Real;
DSTow2Mat = IterateMat(DavidDat.Xtp,DavidDat.Ya,DavidDat.TowMat,2);

TransDataReg = load('./Data/Gen/TransEnergy39Reg');
IndexMatchTp = find([TransDataReg.TpA] == Tp) ; 

DataIn = sqrt(2*TransDataReg.CalTra(IndexMatchTp));
DataOut = sqrt(2*TransDataReg.ExpTra(IndexMatchTp));

figure();
set(gca,'FontSize',30) 
plot(DavidDat.Ya,DST2Mat.*DavidDat.Ya,'-r','DisplayName','WikiWaves - Tranmission','LineWidth',2)
hold on;
plot(DavidDat.Ya,DSTow2Mat.*DavidDat.Ya,'-k','DisplayName','WikiWaves - Tranmission + Overwash' ,'LineWidth',2)

plot(LukeDat.Ya,LBT2Mat.*LukeDat.Ya,'-b','DisplayName','Lukes Code - Tranmission','LineWidth',2 )
plot(LukeDat.Ya,LBTow2Mat.*LukeDat.Ya,'-g','DisplayName','Lukes Code - Tranmission + Overwash','LineWidth',2 )
plot(DataIn,DataOut,'.k','DisplayName', 'Regular Data','MarkerSize',20)
xlabel('amplitude in (a)')
ylabel('amplitude out (a)')
legend()
title(['Comparison Of Linear Models Tp = 0.95'])



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


