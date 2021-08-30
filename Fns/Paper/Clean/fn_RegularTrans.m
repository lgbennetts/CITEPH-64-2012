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

function fn_RegularTransComparison()
% if ~exist('conc','var') conc =79; end %39,79 or empty (1)
% if ~exist('TestNames','var') TestNames ={'6','7','8','9','10','11','11a','11b','12','13','14','15','16','17','18'} ; end 

if ~exist('conc','var') conc =39; end %39,79 or empty (1)
% if ~exist('TestNames','var') TestNames ={'7'}; end 
if ~exist('TestNames','var') TestNames ={'1', '2', '2a', '2b', '3','4', '5', '6', '7'}; end 

%'3'
%  if ~exist('TestNames','var') TestNames ={ '2'}; end 

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

% ExpRef
% ExpTra
% CalRef
% CalTra
proberef = 1:10;
probetra = 11:20;
% OutputStr = 'TargetSpectra Mean_AllProbe_Spectra';
OutputStr = 'TargetSpectra';
for j = 1 : length(TestNames)    
    [CalibTestName,CalFac,Tp,Hs,WaveType,Success]  = fn_FindCalibration(conc,TestNames{j});
%     Success
    outRef = fn_SpectraAverageTgT(Tp,conc,TestNames{j},WaveType,proberef,OutputStr);
    IndexSatisfy =  find(ismember({outRef.name},'SpecFilt')) ;
    ExpRef(j) = outRef(IndexSatisfy).value;  
    
    outTrans = fn_SpectraAverageTgT(Tp,conc,TestNames{j},WaveType,probetra,OutputStr);
    IndexSatisfy =  find(ismember({outTrans.name},'SpecFilt')) ;
    ExpTra(j) = outTrans(IndexSatisfy).value;  

    %Calibration
    outRefC = fn_SpectraAverageTgT(Tp,0,CalibTestName,WaveType,proberef,OutputStr);
    IndexSatisfy =  find(ismember({outTrans.name},'SpecFilt')) ;
    CalRef(j) = (CalFac)^2*outRefC(IndexSatisfy).value; 
    
    outTransC = fn_SpectraAverageTgT(Tp,0,CalibTestName,WaveType,probetra,OutputStr);
    IndexSatisfy =  find(ismember({outTransC.name},'SpecFilt')) ;
    CalTra(j) = (CalFac)^2*outTransC(IndexSatisfy).value; 

    Fac(j) = (CalFac)^2;

    TransConcCal(j) = ExpTra(j) /CalTra(j);
    TransRawConc(j) = ExpTra(j) / ExpRef(j);
    TransCorConc(j) = (2*TransRawConc(j)) / (1 + TransRawConc(j));

    TpA(j) = Tp;
    HsA(j) = Hs;
  
end

MatFile_NM = strcat('Data/Gen/LongNewTransEnergy',int2str(conc),WaveType(1:3));
save(MatFile_NM,'TpA','HsA','Fac', 'ExpRef','ExpTra','CalRef','CalTra','TransConcCal','TransRawConc','TransCorConc');


%figure individual probes spectra - Transmission probes looking weird
% 
% TransSpecMat = outTrans(1).value;
% figure();
% hold on;
% set(gca,'FontSize',18) 
% PlotNames = {};
% for i = 1: size(TransSpecMat,2)
%     
%     plot(outTrans(3).value,TransSpecMat(:,i))
%     PlotNames{end +1} = ['Probe ',num2str(10 + i)];
% end
% 
% plot(outRef(3).value, outRef(2).value , '-r','LineWidth',2);
% PlotNames{end +1} = 'Reflection Probes Disks';
% plot(outTrans(3).value, outTrans(2).value , '-.r','LineWidth',2);
% PlotNames{end +1} = 'Transmission Probes Disks';
% 
% plot(outRefC(3).value, outRefC(2).value , '-b','LineWidth',2);
% PlotNames{end +1} = 'Reflection Probes No Disks';
% plot(outTransC(3).value, outTransC(2).value , '-.b','LineWidth',2);
% PlotNames{end +1} = 'Transmission Probes No Disks';
% 
% xlabel('Period');
% ylabel('Spectra');
% title(['Spectra']);
% legend(PlotNames);

%figure of spectra - lots of transmission?
% figure();
% hold on;
% set(gca,'FontSize',18) 
% clear PlotNames;
% PlotNames = {};
% 
% plot(outRef(2).value, outRef(1).value , '-r');
% PlotNames{end +1} = 'Reflection Probes Disks';
% plot(outTrans(2).value, outTrans(1).value , '-.r');
% PlotNames{end +1} = 'Transmission Probes Disks';
% 
% plot(outRefC(2).value, outRefC(1).value , '-b');
% PlotNames{end +1} = 'Reflection Probes No Disks';
% plot(outTransC(2).value, outTransC(1).value , '-.b');
% PlotNames{end +1} = 'Transmission Probes No Disks';
 
% xlabel('Period');
% ylabel('Spectra');
% title(['Spectra']);
% legend(PlotNames);

% if conc==39
% dum_c = 100*pi*(0.495^2)/2;
% elseif conc==79
% dum_c = 100*pi*(0.495^2);
% end
% 
% %2Demm prediction
% TwoDEMM = Main_AttnModels(model_pers,dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
% TwoDEMM = (TwoDEMM.value);
% 
% Description = ['Concentration ' num2str(conc) '% '];
% 
% figure();
% hold on;
%  set(gca,'FontSize',18) 
% clear PlotNames;
% PlotNames = {};
% 
% plot(TpA(Fac == 1), TransConcCal(Fac == 1) , '.r' ,'MarkerSize', 16);
% pname = ['Regular Tradition T_a Existing Calibration' ];
% PlotNames{end +1} = pname;
% 
% if ~isempty(isempty(find(Fac ~= 1,1)))
% 
%     plot(TpA(Fac ~= 1), TransConcCal(Fac ~= 1) , 'or' ,'MarkerSize', 16);
%     pname = ['Regular Tradition T_a Scaled Other Calibration' ];
%     PlotNames{end +1} = pname;
% 
% end
% 
% plot(TpA, TransRawConc , 'xb' ,'MarkerSize', 16);
% pname = ['Raw Trans/Refs T_a' ];
% PlotNames{end +1} = pname;
% 
% plot(TpA, TransCorConc , '^k' ,'MarkerSize', 16);
% pname = ['Corrected Trans/Refs T_a' ];
% PlotNames{end +1} = pname;
%     
% plot(model_pers, TwoDEMM, '--k' , 'LineWidth', 3);
% PlotNames{end +1} = '2dEMM';
% 
% xlabel('Period');
% ylabel('Transmission Coefficient ');
% title(['Transmission Coefficient ', Description]);
% legend(PlotNames);
% xlim([0 , 2 ]);
% ylim([0 , 1 ]);

%Save Matrices
% MatFile_NM = strcat('Data/Gen/NewATrans',int2str(conc),WaveType(1:3));
% save(MatFile_NM,'TpA','HsA','Fac', 'ExpRef','ExpTra','CalRef','CalTra','TransConcCal','TransRawConc','TransCorConc');


return



