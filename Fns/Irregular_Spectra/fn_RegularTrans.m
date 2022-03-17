% function fn_RegularTrans
%
% DESCRIPTION: Generate the transmissions matrices
%
% INPUTS:
%
% OUTPUTS:
%     MatFile_NM - saved as a matrix file
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by

function fn_RegularTrans()
if ~exist('conc','var') conc =79; end %39,79 or empty (1)
% if ~exist('TestNames','var') TestNames ={'6','7','8','9','10','11','11a','11b','12','13','14','15','16','17','18'} ; end 

% if ~exist('conc','var') conc =39; end %39,79 or empty (1)
% if ~exist('TestNames','var') TestNames ={'7'}; end 
% if ~exist('TestNames','var') TestNames ={'1', '2', '2a', '2b', '3','4', '5', '6', '7'}; end 
% if ~exist('TestNames','var') TestNames ={'4', '5', '6', '7'}; end 

%'3'
if conc == 39
    TestNames ={'6','7','8','9','10','11','11a','11b','12','13','14','15','16','17','18'};
%     TestName= {{'14','9'},{'15','10'},{'16','8'}};
elseif conc == 79
    TestNames ={'1', '2', '2a', '2b', '3','4', '5', '6', '7'};
end

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
%     IndexSatisfy =  find(ismember({outRef.name},'SpecFilt')) ;
    ExpRef_A(j) = outRef(2).value;  
    ExpRef_S(j) = outRef(3).value;  
    
    outTrans = fn_SpectraAverageTgT(Tp,conc,TestNames{j},WaveType,probetra,OutputStr);
%     IndexSatisfy =  find(ismember({outTrans.name},'SpecFilt')) ;
    ExpTra_A(j) = outTrans(2).value;  
    ExpTra_S(j) = outTrans(3).value;   

    %Calibration
    outRefC = fn_SpectraAverageTgT(Tp,0,CalibTestName,WaveType,proberef,OutputStr);
    IndexSatisfy =  find(ismember({outTrans.name},'SpecFilt')) ;
%     CalRef(j) = (CalFac)^2*outRefC(IndexSatisfy).value; 
%     CalRef_A(j) = (CalFac)^2*outRefC(2).value;  
%     CalRef_S(j) = (CalFac)^2*outRefC(3).value;  
    CalRef_A(j) = (CalFac)*outRefC(2).value;  
    CalRef_S(j) = (CalFac)*outRefC(3).value;  
    
    outTransC = fn_SpectraAverageTgT(Tp,0,CalibTestName,WaveType,probetra,OutputStr);
    IndexSatisfy =  find(ismember({outTransC.name},'SpecFilt')) ;
%     CalTra(j) = (CalFac)^2*outTransC(IndexSatisfy).value;
    CalTra_A(j) = (CalFac)*outTransC(2).value;  
    CalTra_S(j) = (CalFac)*outTransC(3).value;  

    Fac(j) = (CalFac);

    TransConcCal_A(j) = ExpTra_A(j) /CalTra_A(j);
    TransRawConc_A(j) = ExpTra_A(j) / ExpRef_A(j);
    TransCorConc_A(j) = (2*TransRawConc_A(j)) / (1 + TransRawConc_A(j));

    TransConcCal_S(j) = ExpTra_S(j) /CalTra_A(j);
    TransRawConc_S(j) = ExpTra_S(j) / ExpRef_A(j);
    TransCorConc_S(j) = (2*TransRawConc_S(j)) / (1 + TransRawConc_A(j));    
    
    TpA(j) = Tp;
    HsA(j) = Hs;
  
end

MatFile_NM = strcat('Data/Gen/NewTransAmp',int2str(conc),WaveType(1:3));
save(MatFile_NM,'TpA','HsA','Fac', 'ExpRef_A','ExpTra_A','CalRef_A','CalTra_A','TransConcCal_A','TransRawConc_A','TransCorConc_A', 'ExpRef_S','ExpTra_S','CalRef_S','CalTra_S','TransConcCal_S','TransRawConc_S','TransCorConc_S');


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



