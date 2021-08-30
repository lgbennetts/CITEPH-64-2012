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

function fn_TransmissionReg(conc,TestName,probes)

if ~exist('PLOT_TYP','var'); PLOT_TYP ='trans_coeff'; end

if ~exist('conc','var'); conc= 79; end ;

if ~exist('binnum','var');  binnum=30; end;

if ~exist('Cols','var'); Cols= '#058000'; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;


%Load data
if conc == 39
    TransData = load('./Data/Trans38');
else
    TransData = load('./Data/Trans77');
end

Alpha = -log(TransData.trans)./5;

Description = ['Concentration ' num2str(conc) '% '];

% Models
if conc==39
dum_c = 100*pi*(0.495^2)/2;
elseif conc==79
dum_c = 100*pi*(0.495^2);
end

%Predictions
Boltzman0 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
Boltzman0 = (Boltzman0.value).^2;
Boltzman0 = Boltzman0.';
Boltzman0 = [0;Boltzman0;1];

Boltzman1 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
Boltzman1 = (Boltzman1.value).^2;
Boltzman1 = Boltzman1.';
Boltzman1 = [0;Boltzman1;1];

%2Demm prediction
TwoDEMM = Main_AttnModels(model_pers(2:end-1),dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
TwoDEMM = (TwoDEMM.value).^2;
TwoDEMM = TwoDEMM.';
TwoDEMM = [0;TwoDEMM;1];

ffmod = 1.0./ model_pers.';
ffmod = ffmod(2:end-1);
TEM = TwoDEMM(2:end-1);
TB0 = Boltzman0(2:end-1);
TB1 = Boltzman1(2:end-1);
AEM = -(log(TEM)./5);
AB1 = -(log(TB1)./5);
AB0 = -(log(TB0)./5);

%Transmission PLots
%Spectra Comparisons - both Calibration and Experimental


%Attenuation Coefficient
%What if we just fit one?
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
fflm  = 0.1:0.1:10;
PlotNames = {};

plot(1./TransData.T, TransData.trans , '.','Color',Cols ,'MarkerSize', 16);
pname = ['Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

Upper = TransData.trans_std(1,:);
Lower = TransData.trans_std(2,:);
ffrem0 = 1./TransData.T;
ShadeX = [ffrem0; flipud(ffrem0)];
inBetween = [Upper ; flipud(Lower)];
patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols)
alpha(0.3)
PlotNames{end +1} = 'Max and Min Region (Among Probes)';

plot(ffrem0, Upper, ':','Color',Cols , 'LineWidth', 1);
PlotNames{end +1} = 'Max and Min Boundary';
plot(ffrem0, Lower, ':','Color',Cols , 'LineWidth', 1);
PlotNames{end +1} = '';
    

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


figure();
plot(1./TransData.T, Alpha , '.','Color',Cols,'MarkerSize', 16);
hold on;
PlotNames = {};
pname = ['Conc ', num2str(conc),' |  Transmission Probes '];
PlotNames{end +1} = pname;

y = log(Alpha(Alpha > 0));
x = log(1./TransData.T(Alpha > 0)); 


mdl = fitlm(x,y);  
b = mdl.Coefficients.Estimate(1);
m = mdl.Coefficients.Estimate(2); 
plot(exp(x), exp(m*x + b),':','Color',Cols , 'LineWidth', 2);
pname = ['Model Slope : ', num2str(m)];
PlotNames{end +1} = pname;

plot(ffmod, AB0, '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(ffmod, AB1, '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(ffmod, AEM, '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

xlim([0.5 , 2.55 ]);

xlabel('Frequency');
ylabel('Attenuation Coefficient');
title(['Attenuation Coefficient', Description]);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend(PlotNames);


% title(['Transmission Coefficients Between 3Tp/4 and 6Tp/4 ' ,Description]);
% xlim([0.4,2]);
% ylim([0,1]);
% xlabel('Period [s]')
% ylabel('Transmission Coefficient')
% legend();

% figure();
% plot(1./ff1{1},MovTransmissionsCoeff{1},'xk','DisplayName',Lab, 'MarkerSize',5);
% hold on;
% plot(model_pers,Boltzman0,'--r','DisplayName','Boltzmann Remove Scattering','LineWidth',2);
% plot(model_pers,Boltzman1,'--b','DisplayName','Boltzmann','LineWidth',2);
% 
% plot(model_pers,TwoDEMM,'-.k','DisplayName','2D EMM','LineWidth',2);
% 
% title(['Moving Average (7) Transmission Coefficients - Ratio of Binned Averages of Spectral Energy, number of bins : ' num2str(binnum)]);
% xlim([0.4,2]);
% ylim([0,1]);
% legend();

return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

