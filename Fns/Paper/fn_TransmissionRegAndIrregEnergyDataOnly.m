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

if ~exist('ColsReg','var'); ColsReg= 'green'; end ;
if ~exist('Symbols','var'); Symbols= {'x','d','^','s'}; end ;
if ~exist('Cols','var');Cols= {'blue','red','black','cyan'}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=2e2; end
if ~exist('model_pers','var'); model_pers=0.5:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=1; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;


%Load data
if conc == 39
%     TransDataReg = load('./Data/Trans38');
    TransDataReg = load('./Data/Gen/TransEnergy39Reg');
    TransDataIrreg = load('./Data/Gen/TransEnergy39Irr');
else
%     TransDataReg = load('./Data/Trans77');
    TransDataReg = load('./Data/Gen/TransEnergy79Reg');
    TransDataIrreg = load('./Data/Gen/TransEnergy79Irr');
end

%sqrt for amplitude, and calculate wl
for i = 1:length(TransDataIrreg.Trans)
%     TransDataIrreg.Trans{i} = sqrt(TransDataIrreg.Trans{i});
    
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





Description = ['Concentration ' num2str(conc) '% '];



%Wavelength
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};


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
ylabel('Transmission Coefficient (Energy) ');
title(['Transmission Coefficient (Energy) ', Description]);
legend(PlotNames);
xlim([0 , 5 ]);
ylim([0 , 1 ]);

%3D plot
Size = 60;
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};

XData = WLReg;
YData = sqrt(2*TransDataReg.CalRef);
ZData = TransDataReg.TransConcCal;

scatter3(WLReg,sqrt(2*TransDataReg.CalRef),TransDataReg.TransConcCal,Size,'green','filled')
pname = ['Regular Conc ', num2str(conc),' |  Tranmission Probes'];
PlotNames{end +1} = pname;

for j = 1: length(TransDataIrreg.TpTarg)
    XData = [XData,TransDataIrreg.WL{j}.'];
    YData = [YData,sqrt(2*TransDataIrreg.Ain{j}).'];
    ZData = [ZData,TransDataIrreg.Trans{j}.'];

    scatter3(TransDataIrreg.WL{j},sqrt(2*TransDataIrreg.Ain{j}), TransDataIrreg.Trans{j} ,Size,Cols{j} ,'filled');
    pname = ['Irregular Conc ', num2str(conc),' |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end
    

xlabel('Wavelength (m)');
ylabel('Incoming Amplitude (m)');
zlabel('Transmission Coefficient (Energy) ');


%3D function
[xq,yq] = meshgrid(0:0.1:5.4,0:0.01:0.07);
z1 = griddata(XData,YData,ZData,xq,yq,'v4');
% plot3(XData,YData,ZData,'mo')
surf(xq,yq,z1)
PlotNames{end +1} = 'Biharmonic Splines';  

xlabel('Wavelength (m)');
ylabel('Incoming Amplitude (m)');
zlabel('Transmission Coefficient (Energy) ');
legend(PlotNames);



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

