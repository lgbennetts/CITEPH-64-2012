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

if ~exist('conc','var'); conc= 39; end ; %79

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
if conc == 39
%     TransDataReg = load('./Data/Trans38');
    TransDataReg = load('./Data/Gen/NewATrans39Reg');
    TransDataIrreg = load('./Data/Gen/NewTrans39Irr');
else
%     TransDataReg = load('./Data/Trans77');
    TransDataReg = load('./Data/Gen/NewATrans79Reg');
    TransDataIrreg = load('./Data/Gen/NewTrans79Irr');
end

%sqrt for amplitude, and calculate wl
for i = 1:length(TransDataIrreg.Trans)
    TransDataIrreg.Trans{i} = sqrt(TransDataIrreg.Trans{i});
    
    for j = 1:length(TransDataIrreg.WL{i})
        ki = dispersion_free_surface((2*pi*TransDataIrreg.F{i}(j))^2/9.81,0,3.1);
        TransDataIrreg.WLN{i}(j,1) = 2*pi/imag(ki);
    end
    
end

%Reg
for i = 1: length(TransDataReg.TpA)
    ki = dispersion_free_surface((2*pi/TransDataReg.TpA(i))^2/9.81,0,3.1);
    WLReg(i) = 2*pi/imag(ki);
end


%Alphas
AlphaReg = -(log(TransDataReg.TransConcCal)./5);
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
TwoDEMM = TwoDEMM.value;%(TwoDEMM.value).^2;
TwoDEMM = TwoDEMM.';
TwoDEMM = [0;TwoDEMM;1];


%Predictions
Boltzman0 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
% Boltzman0 = (Boltzman0.value).^2;
Boltzman0 = (Boltzman0.value);
Boltzman0 = Boltzman0.';
Boltzman0 = [0;Boltzman0;1];

Boltzman1 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
% Boltzman1 = (Boltzman1.value).^2;
Boltzman1 = (Boltzman1.value);
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
PlotNames{end +1} = '2dEMM';

%Regular
plot(WLReg, TransDataReg.TransConcCal , '.','Color',ColsReg ,'MarkerSize', 25);
pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

%Irregular
for j = 1: length(TransDataIrreg.TpTarg)
    plot(TransDataIrreg.WL{j}, TransDataIrreg.Trans{j} , Symbols{j},'Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end
    


xlabel('Wavelength (m)');
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

%Models
plot(wlMod.', AB0.', '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(wlMod.', AB1.', '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(wlMod.', AEM.', '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

%Regular
plot(WLReg, AlphaReg , '.','Color',ColsReg ,'MarkerSize', 25);
pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

xl = [];
yl = [];

%LR
y = log(AlphaReg(AlphaReg > 0));
x = log(WLReg(AlphaReg > 0)); 

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
    plot(TransDataIrreg.WL{j}, AlphaIrreg{j} , Symbols{j},'Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
    
    %LR
    y = log(AlphaIrreg{j}(AlphaIrreg{j} > 0));
    x = log(TransDataIrreg.WL{j}(AlphaIrreg{j} > 0)); 

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


xlabel('Wavelength (m)');
ylabel('Attenuation Coefficient');
title(['Attenuation Coefficient', Description]);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend(PlotNames);



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

