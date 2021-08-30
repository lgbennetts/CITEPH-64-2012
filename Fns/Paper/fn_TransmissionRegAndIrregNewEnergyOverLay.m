% function fn_InvestigateSpectra
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
%
% General:
%
% conc      = either 39 or 79 , otherwise calibration
% WaveType =  'Regular' or 'Irregular'
% TestName = name of test in fn_*_testspecs function for asociated
% experiment
% probe - probe to read, allow any from  1- 20 (10 on left, and 10 on right of probe)
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function fn_TransmissionRegAndIrreg(conc,TestName,probes)

if ~exist('PLOT_TYP','var'); PLOT_TYP ='trans_coeff'; end

if ~exist('conc','var'); conc= 79; end ; %79

if ~exist('binnum','var');  binnum=30; end;

if ~exist('ColsReg','var'); ColsReg= '#058000'; end ;
if ~exist('Symbols','var'); Symbols= {'x','d','^','s'}; end ;
if ~exist('Cols','var');Cols= {'#ff0000','#0bff01 ','#0487f9','#9701ff'}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.5:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;


%Load data
TransDataReg39 = load('./Data/Gen/TransEnergy39Reg');
TransDataIrreg39 = load('./Data/Gen/TransEnergy39Irr');
TransDataReg79 = load('./Data/Gen/TransEnergy79Reg');
TransDataIrreg79 = load('./Data/Gen/TransEnergy79Irr');

%sqrt for amplitude, and calculate wl
for i = 1:length(TransDataIrreg39.Trans)
%     TransDataIrreg.Trans{i} = sqrt(TransDataIrreg.Trans{i});
    
    for j = 1:length(TransDataIrreg39.WL{i})
        ki = dispersion_free_surface((2*pi*TransDataIrreg39.F{i}(j))^2/9.81,0,3.1);
        TransDataIrreg39.WLN{i}(j,1) = 2*pi/imag(ki);
    end
    
end

for i = 1:length(TransDataIrreg79.Trans)
%     TransDataIrreg.Trans{i} = sqrt(TransDataIrreg.Trans{i});
    
    for j = 1:length(TransDataIrreg79.WL{i})
        ki = dispersion_free_surface((2*pi*TransDataIrreg79.F{i}(j))^2/9.81,0,3.1);
        TransDataIrreg79.WLN{i}(j,1) = 2*pi/imag(ki);
    end
    
end

%Reg
for i = 1: length(TransDataReg39.TpA)
    ki = dispersion_free_surface((2*pi/TransDataReg39.TpA(i))^2/9.81,0,3.1);
    WLReg39(i) = 2*pi/imag(ki);
end

for i = 1: length(TransDataReg79.TpA)
    ki = dispersion_free_surface((2*pi/TransDataReg79.TpA(i))^2/9.81,0,3.1);
    WLReg79(i) = 2*pi/imag(ki);
end


%Alphas
% AlphaReg = -(log(TransDataReg.TransConcCal)./5);
% for j = 1: length(TransDataIrreg.TpTarg)
%     AlphaIrreg{j}  = -log(TransDataIrreg.Trans{j})/5;
% end

Description = ['Both Concentrations% '];

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
    wlMod(i) = 2*pi/(imag(ki));
end



%Wavelength
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
fflm  = 0.1:0.1:10;
PlotNames = {};

plot(wlMod.', TB0.', '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(wlMod.', TB1.', '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(wlMod.', TEM.', '--k' , 'LineWidth', 3);
PlotNames{end +1} = 'Lukes';

%Regular
plot(WLReg39, (TransDataReg39.TransConcCal) , '.','Color',ColsReg ,'MarkerSize', 25);
pname = ['Regular Conc 39', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

%Irregular
for j = 1: length(TransDataIrreg39.TpTarg)
    plot(TransDataIrreg39.WL{j}, sqrt(TransDataIrreg39.Trans{j}) , Symbols{j},'Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc 39 ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg39.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg39.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end

plot(WLReg79, (TransDataReg79.TransConcCal) , '.','Color',ColsReg ,'MarkerSize', 25);
pname = ['Regular Conc 79', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

%Irregular
for j = 1: length(TransDataIrreg79.TpTarg)
    plot(TransDataIrreg79.WL{j}, sqrt(TransDataIrreg79.Trans{j}) , Symbols{j},'Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc 79 ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg79.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg79.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end
    


xlabel('Wavelength (m)');
ylabel('Transmission Coefficient (Energy) ');
title(['Transmission Coefficient (Energy) ', Description]);
legend(PlotNames);
xlim([0 , 5 ]);
ylim([0 , 1 ]);



% %Attenuation Coefficient
% %What if we just fit one?
% figure();
% hold on;
% set(gca,'FontSize',18) 
% clear PlotNames;
% PlotNames = {};
% 
% %Models
% plot(wlMod.', AB0.', '-.r' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Conservative';
% 
% plot(wlMod.', AB1.', '-.b' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Losing Energy';
% 
% plot(wlMod.', AEM.', '--k' , 'LineWidth', 3);
% PlotNames{end +1} = 'Lukes';
% 
% %Regular
% plot(WLReg, AlphaReg , '.','Color',ColsReg ,'MarkerSize', 25);
% pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
% PlotNames{end +1} = pname;
% 
% xl = [];
% yl = [];
% 
% %LR
% y = log(AlphaReg(AlphaReg > 0));
% x = log(WLReg(AlphaReg > 0)); 
% 
% xl = [xl;x.'];
% yl = [yl;y.'];
% 
% mdl = fitlm(x,y);  
% b = mdl.Coefficients.Estimate(1);
% m = mdl.Coefficients.Estimate(2); 
% plot(exp(x), exp(m*x + b),':','Color',ColsReg , 'LineWidth', 2);
% pname = ['LinearRegress Slope : ', num2str(m)];
% PlotNames{end +1} = pname;
% %Irregular
% for j = 1: length(TransDataIrreg.TpTarg)
%     plot(TransDataIrreg.WL{j}, AlphaIrreg{j} , Symbols{j},'Color',Cols{j} ,'MarkerSize', 16);
%     pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
%     PlotNames{end +1} = pname;  
%     
%     %LR
%     y = log(AlphaIrreg{j}(AlphaIrreg{j} > 0));
%     x = log(TransDataIrreg.WL{j}(AlphaIrreg{j} > 0)); 
% 
%     xl = [xl;x];
%     yl = [yl;y];
% 
%     mdl = fitlm(x,y);  
%     b = mdl.Coefficients.Estimate(1);
%     m = mdl.Coefficients.Estimate(2); 
%     plot(exp(x), exp(m*x + b),':','Color',Cols{j} , 'LineWidth', 2);
%     pname = ['LinearRegress Slope : ', num2str(m)];
%     PlotNames{end +1} = pname;
% end
%     
% mdl = fitlm(xl,yl);  
% b = mdl.Coefficients.Estimate(1);
% m = mdl.Coefficients.Estimate(2); 
% plot(exp(xl), exp(m*xl + b),':k', 'LineWidth', 2);
% pname = ['LinearRegress Slope : ', num2str(m)];
% PlotNames{end +1} = pname;
% 
% 
% xlabel('Wavelength (m)');
% ylabel('Attenuation Coefficient');
% title(['Attenuation Coefficient', Description]);
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend(PlotNames);



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

