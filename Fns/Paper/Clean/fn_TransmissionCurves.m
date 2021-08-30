
close all;


TransDataReg = load('./Data/Gen/NewTransEnergy39Reg.mat');
TransDataIrreg = load('./Data/Gen/NewTransEnergy39Irr.mat');


%Reg
for i = 1: length(TransDataReg.TpA)
    ki = dispersion_free_surface((2*pi/TransDataReg.TpA(i))^2/9.81,0,3.1);
    WLReg(i) = 2*pi/imag(ki);
end


%Alphas
AlphaReg = -(log(TransDataReg.trans)./5);
for j = 1: length(TransDataIrreg.TpTarg)
    AlphaIrreg{j}  = -log(TransDataIrreg.Trans{j})/5;
end

Description = ['Concentration ' num2str(conc) '% '];

% Models
if conc==39
dum_c = 100*pi*(0.495^2)/2;
elseif conc==79
dum_c = 100*pi*(0.495^2);
end

%2Demm prediction
TwoDEMM = Main_AttnModels(model_pers(2:end-1),dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
TwoDEMM = (TwoDEMM.value).^2;
TwoDEMM = TwoDEMM.';
TwoDEMM = [0;TwoDEMM;1];


%Predictions
Boltzman0 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
Boltzman0 = (Boltzman0.value).^2;
Boltzman0 = Boltzman0.';
Boltzman0 = [0;Boltzman0;1];

Boltzman1 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
Boltzman1 = (Boltzman1.value).^2;
Boltzman1 = Boltzman1.';
Boltzman1 = [0;Boltzman1;1];

ffmod = 1.0./ model_pers.';
ffmod = ffmod(2:end-1);
TEM = TwoDEMM(2:end-1);
TB0 = Boltzman0(2:end-1);
TB1 = Boltzman1(2:end-1);
AEM = -(log(TEM)./5);
AB1 = -(log(TB1)./5);
AB0 = -(log(TB0)./5);

for i = 1: length(ffmod)
    ki = dispersion_free_surface((2*pi*ffmod(i))^2/9.81,0,3.1);
    wlMod(i) = 2*pi/imag(ki);
end


%Transmission PLots
%Spectra Comparisons - both Calibration and Experimental
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
fflm  = 0.1:0.1:10;
PlotNames = {};

%Regular
plot(1./TransDataReg.T, TransDataReg.trans , '.','Color',ColsReg ,'MarkerSize', 16);
pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

%Irregular
for j = 1: length(TransDataIrreg.TpTarg)
    plot(TransDataIrreg.F{j}, TransDataIrreg.Trans{j} , '.','Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end
    

plot(ffmod.', TB0.', '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(ffmod.', TB1.', '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(ffmod.', TEM.', '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

xlabel('Frequency');
ylabel('Transmission Coefficient ');
title(['Transmission Coefficient ', Description]);
legend(PlotNames);
xlim([0 , 5 ]);
ylim([0 , 1 ]);

%wavelength
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
fflm  = 0.1:0.1:10;
PlotNames = {};

%Regular
plot(WLReg, TransDataReg.trans , '.','Color',ColsReg ,'MarkerSize', 16);
pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

%Irregular
for j = 1: length(TransDataIrreg.TpTarg)
    plot(TransDataIrreg.WL{j}, TransDataIrreg.Trans{j} , '.','Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end
    

plot(wlMod.', TB0.', '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(wlMod.', TB1.', '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(wlMod.', TEM.', '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

xlabel('Wave Length');
ylabel('Transmission Coefficient ');
title(['Transmission Coefficient ', Description]);
legend(PlotNames);
xlim([0 , 5 ]);
ylim([0 , 1 ]);



%Attenuation Coefficient
%What if we just fit one?
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};

%Regular
plot(1./TransDataReg.T, AlphaReg , '.','Color',ColsReg ,'MarkerSize', 16);
pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

xl = [];
yl = [];

%LR
y = log(AlphaReg(AlphaReg > 0));
x = log(1./TransDataReg.T(AlphaReg > 0)); 

xl = [xl;x.'];
yl = [yl;y.'];

mdl = fitlm(x,y);  
b = mdl.Coefficients.Estimate(1);
m = mdl.Coefficients.Estimate(2); 
plot(exp(x), exp(m*x + b),':','Color',ColsReg , 'LineWidth', 2);
pname = ['LinearRegress Slope : ', num2str(m)];
PlotNames{end +1} = pname;
%Irregular
for j = 1: length(TransDataIrreg.TpTarg)
    plot(TransDataIrreg.F{j}, AlphaIrreg{j} , '.','Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
    
    %LR
    y = log(AlphaIrreg{j}(AlphaIrreg{j} > 0));
    x = log(TransDataIrreg.F{j}(AlphaIrreg{j} > 0)); 

    xl = [xl;x];
    yl = [yl;y];

    mdl = fitlm(x,y);  
    b = mdl.Coefficients.Estimate(1);
    m = mdl.Coefficients.Estimate(2); 
    plot(exp(x), exp(m*x + b),':','Color',Cols{j} , 'LineWidth', 2);
    pname = ['LinearRegress Slope : ', num2str(m)];
    PlotNames{end +1} = pname;
end
    
mdl = fitlm(xl,yl);  
b = mdl.Coefficients.Estimate(1);
m = mdl.Coefficients.Estimate(2); 
plot(exp(xl), exp(m*xl + b),':k', 'LineWidth', 2);
pname = ['LinearRegress Slope : ', num2str(m)];
PlotNames{end +1} = pname;

plot(ffmod.', AB0.', '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(ffmod.', AB1.', '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(ffmod.', AEM.', '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

xlabel('Frequency');
ylabel('Attenuation Coefficient');
title(['Attenuation Coefficient', Description]);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend(PlotNames);



