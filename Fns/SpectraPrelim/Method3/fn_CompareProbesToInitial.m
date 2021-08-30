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

function fn_CompareProbesToInitial(conc, probes,TestName)


%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Which user?
[user_number,user_name] = citeph_get_user();

if user_number == 0
    OTHER_USR=1;
end


%% GENERAL
%pick expriment

%Have to change the way we pull data, to get same Tp, Hs

if ~exist('conc','var') conc=39; end %39,79 or empty (1)
if ~exist('probe','var');     probes= 11:20 ;end; 
if ~exist('TestName','var');     TestName= '9' ;end; 
if ~exist('TargetRegion','var');  TargetRegion=[0.6,1.2]; end



% if ~exist('data_out','var')
%  data_out.name={'amp-harmo-steady-1st'};
%  data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
% end
% 
% if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

%Read Data in Probe 1

%Find test with concentration,Tp,Hs
%As well as if it has a calibration

TestNames = {'1','2', '2a','2b','3','4','5','6','7','8','9','10'};
%TestNames = {'9'};

t0description ='waves reach x';
t1description ='final waves reach x';
%t1description ='beach ref waves reach x';
    
    
% TestNames = {'1'};
trans_std = [];
T = [];
trans = [];
H = [];
Types = [];
Index = [];
for TestName = TestNames

    [test_pram,test_WaveType,test_Success,calib_Success,calib_Name] ...
                = fn_FindTestAndCalib(conc,TestName);
    if (test_Success == 1 && calib_Success ==1 )

        calibconc = 1;
        allprobes = 1:20;
        Tp = test_pram.period/10.0;
        Hs = 10*test_pram.wave_height;
        
        TargetRegion=[0.55*Tp,1.95*Tp];

        [Tps,TpsRegions] = fn_GetPeaksIncoming(calibconc ,calib_Name,1:10,Tp,TargetRegion,t0description,t1description);

        for i = 1:length(Tps)
            TpI = Tps(i);
            TpsRegionI = TpsRegions(i,:);

            AvgEnergyInc = fn_ExtractEnergy(calibconc,calib_Name,probes,TpI,TpsRegionI,conc,t0description,t1description);
            AvgEnergyOut = fn_ExtractEnergy(conc,TestName,probes,TpI,TpsRegionI,conc,t0description,t1description);

            TransMean(i) = AvgEnergyOut/ AvgEnergyInc;
        end
        trans_std = [trans_std,[TransMean;TransMean] ];
        T = [T,Tps];
        trans = [trans,TransMean ];
        H = [H, Hs*ones(1,length(TransMean))];
        
        Type    = cell(1, length(TransMean));
        Type(:) = {test_WaveType};
        Types = [Types,Type ];
        
        Index = [Index,1:length(Tps) ];
        clear TransMean Type;

    end
end

MatFile_NM = strcat('Data/Gen/Trans',int2str(conc),'All');
  save(MatFile_NM,'T','trans_std','trans','H','Types','Index');

% [AvgSn,Period] = fn_ReadandCombineSpectra(conc,TestName,Tp,probes);
% 
% %Are these good peaks?
% %plot against FFT
% 
% figure();
% plot(Period,AvgSn,'-');
% hold on;
% for i = 1: length(Tps)
%     plot([Tps(i),Tps(i)],[0,max(AvgSn)],'--k');
%     txt = ['T = ' num2str(Tps(i))];
%     text(Tps(i),0.5*max(AvgSn),txt)
% end
% title([description '  Found Spectral Peaks (From Probe 1) and Spectrogram for Ideal Time Window For Tp']);
% xlabel('Period (s)');
% xlabel('Spectogram ');
% hold off;


%Other Probe

return


function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return

