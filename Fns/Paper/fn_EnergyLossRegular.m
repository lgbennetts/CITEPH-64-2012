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
clear all;
close all;

TransDataReg = load('./Data/Gen/NewTransEnergy39Reg');
TransDataIrr = load('./Data/Gen/NewTransEnergy39Irr');

EnergyExp = (TransDataReg.ExpRef + TransDataReg.ExpTra - TransDataReg.CalTra) ./ (TransDataReg.CalTra);

for i = 1:length(TransDataIrr.CalibRef)
    
    TransDataIrr.EnergyExp{i} = (TransDataIrr.ExpRef{i} + TransDataIrr.ExpTra{i}  -  TransDataIrr.CalibRef{i} ) ./ TransDataIrr.CalibTra{i};
end
 

figure('DefaultAxesFontSize',18);
PlotNames = {};
scatter(TransDataReg.TpA,EnergyExp,80,TransDataReg.HsA/2,'filled')
pname = ['Regular Experiments (I+R) + T - I`'];
PlotNames{end +1} = pname;
h = colorbar();
ylabel(h, 'amplitude in (m)')
title(['Total In-Plane Energy']);
xlabel('Tp (s)')
ylabel('Energy')
axis([0 2 0 1])
legend(PlotNames);

for i = 1:length(TransDataIrr.CalibRef)
    
    figure('DefaultAxesFontSize',18);
    PlotNames = {};
    scatter(1./TransDataIrr.F{i},TransDataIrr.EnergyExp{i},80,sqrt(2*TransDataIrr.CalibTra{i}),'filled')
    pname = ['Irregular Experiments (I+R) + T - I` - Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})];
    PlotNames{end +1} = pname;
    h = colorbar();
    ylabel(h, 'amplitude in (m)')
    title(['Total In-Plane Energy Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})]);
    xlabel('Tp (s)')
    ylabel('Energy')
    xlim([0.25*TransDataIrr.TpTarg{i} 2*TransDataIrr.TpTarg{i}])
    legend(PlotNames);
    
    figure('DefaultAxesFontSize',18);
    PlotNames = {};
    hold on;
    stem(1./TransDataIrr.F{i},TransDataIrr.CalibRef{i})
    pname = ['Calibration Reflection Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})];
    PlotNames{end +1} = pname;
    stem(1./TransDataIrr.F{i},TransDataIrr.CalibTra{i})
    pname = ['Calibration Transmission Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})];
    PlotNames{end +1} = pname;
    stem(1./TransDataIrr.F{i},TransDataIrr.ExpRef{i})
    pname = ['Experiment Reflection Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})];
    PlotNames{end +1} = pname;
    stem(1./TransDataIrr.F{i},TransDataIrr.ExpTra{i})
    pname = ['Experiment Transmission Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})];
    PlotNames{end +1} = pname;
    
    title(['Spectra Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})]);
    xlabel('Tp (s)')
    ylabel('S(T)')
    legend(PlotNames);

end

figure('DefaultAxesFontSize',18);
PlotNames = {};
scatter(TransDataReg.TpA,TransDataReg.TransConcCal,80,TransDataReg.HsA/2,'filled')
pname = ['Regular Experiments'];
PlotNames{end +1} = pname;
h = colorbar();
ylabel(h, 'amplitude in (m)')
title(['Transmission']);
xlabel('Tp (s)')
ylabel('Transmission')
axis([0 2 0 1])
legend(PlotNames);

for i = 1:length(TransDataIrr.CalibRef)
    
    figure('DefaultAxesFontSize',18);
    PlotNames = {};
    scatter(1./TransDataIrr.F{i},TransDataIrr.Trans{i},80,sqrt(2*TransDataIrr.CalibTra{i}),'filled')
    pname = ['Irregular Experiments  T  | Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})];
    PlotNames{end +1} = pname;
    h = colorbar();
    ylabel(h, 'amplitude in (m)')
    title(['Transmission Tp: ',num2str(TransDataIrr.TpTarg{i}), '  Hs ',num2str(TransDataIrr.HsTarg{i})]);
    xlabel('Tp (s)')
    ylabel('Transmission')
    xlim([0.25*TransDataIrr.TpTarg{i} 2*TransDataIrr.TpTarg{i}])
    legend(PlotNames);
   

end


figure('DefaultAxesFontSize',18);
PlotNames = {};
scatter(TransDataReg.TpA,(TransDataReg.ExpRef - TransDataReg.CalTra) ./ TransDataReg.CalTra,80,TransDataReg.HsA/2,'filled')
pname = ['Regular Experiments'];
PlotNames{end +1} = pname;
h = colorbar();
ylabel(h, 'amplitude in (m)')
title(['In Plane Reflection']);
xlabel('Tp (s)')
ylabel('In Plane Reflection')
axis([0 2 -0.1 1])
legend(PlotNames);

% figure('DefaultAxesFontSize',18);
% scatter(TransDataReg.TpA,EnergyExp,80,TransDataReg.HsA/2,'filled','DisplayName','Regular Experiments (I+R) + T - I`')
% pname = ['Regular Transmission'];
% PlotNames{end +1} = pname;
% h = colorbar();
% ylabel(h, 'amplitude in (m)')
% title(['Total Energy']);
% xlabel('Tp (s)')
% ylabel('Energy')
% axis([0 2 0 1])
% legend();


% for i = 1:length(TransDataIrr.CalibRef)
%     plot( TransDataIrr.TpTarg{i},TransDataIrr.EnergyExp{i},'ok', 'MarkerSize',12, 'DisplayName',['Energy ', num2str(TransDataIrr.TpTarg{i}) ,' ', num2str(TransDataIrr.HsTarg{i}) ]);
%     hold on;
% end

    

figure('DefaultAxesFontSize',18);
plot(TransDataReg.HsA,EnergyExp,'ok', 'MarkerSize',12, 'DisplayName','Energy Loss (I+R) + T  - I`');
hold on;
title(['Total Energy']);
xlabel('Hs (m)')
ylabel('Energy')
axis([0 0.08 0 1])
legend();

% %Load data
% if conc == 39
% %     TransDataReg = load('./Data/Trans38');
%     TransDataReg = load('./Data/Gen/TransEnergy39Reg');
%     TransDataIrreg = load('./Data/Gen/TransEnergy39Irr');
% else
% %     TransDataReg = load('./Data/Trans77');
%     TransDataReg = load('./Data/Gen/TransEnergy79Reg');
%     TransDataIrreg = load('./Data/Gen/TransEnergy79Irr');
% end



