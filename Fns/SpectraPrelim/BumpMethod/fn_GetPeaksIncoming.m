% function fn_JustFourierMovingWindow
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
% tm - time series times in seconds
% data - data at times (such as displacement in meters)
% Tpers - string containing function call, to be evaluated to give Tpers.
%
% General:
%
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function [Tps,TpsRegions] = fn_GetPeaksIncoming(conc,TestName,probes,Tp,TargetRegion,t0description,t1description,TWindow)


    if ~exist('conc','var');  conc=1; end
    if ~exist('TestName','var');  TestName='1'; end
    if ~exist('probes','var');  probes=11:20; end
    if ~exist('TargetRegion','var');  TargetRegion=[0.6,1.2]; end
    if ~exist('EnergyFactorTolerance','var');  EnergyFactorTolerance=10.0; end
    
    PeakLim = 30;
    
    if ~exist('Tp','var')
            [test_pram,test_WaveType,test_Success,calib_Success,calib_Name] ...
                = fn_FindTestAndCalib(conc,TestName);
            Tp = test_pram.period/10.0;
    end
    
    [Sn,periods,SnTimes] = fn_ReadandCombineSpectra(conc,TestName,Tp,probes,t0description,t1description,TWindow);
    
    Sn(or(periods < TargetRegion(1) ,periods > TargetRegion(2) )) = 0;
    
%     figure();
%     plot(periods,Sn,'-b','LineWidth',2);
%     [test_pram,test_WaveType,test_Success,calib_Success,calib_Name] = fn_FindTestAndCalib(conc,TestName);
%     title(['Concentation -  ',num2str(conc),'  ',test_WaveType,' Wave ', 'Tp - ',num2str(Tp), '[s] Hs - ',num2str(test_pram.wave_height*10),'[mm]  [' t0description,'(' num2str(SnTimes(1)),'s)' ,' , ',t1description,'(' num2str(SnTimes(2)),'s)' , ']']);
%     xlim([0,2*Tp]);
%     hold on;
    
    [TpC , TpRegion,PeakEnergy ] = fn_ExtractFirstPeak(Sn,periods);
    
    i= 1 ;
    Tps(i) = TpC;
    TpsRegions(i,:) = TpRegion ;
    Energy = PeakEnergy;
    while and(Energy > PeakEnergy/EnergyFactorTolerance, i < PeakLim )
        i = i +1;
        [TpC , TpRegion,Energy ] = fn_ExtractOtherPeaks(Sn,periods,TpsRegions);
        if Energy > PeakEnergy/EnergyFactorTolerance
            Tps(i) = TpC;
            TpsRegions(i,:) = TpRegion ;
        end
    end
    


return

function [TpC,TpRegion,Energy] = fn_ExtractOtherPeaks(An,periods,RemovedRegions)

Sn = An;

%Zero out remove regions
if ~isempty(RemovedRegions)
    for j = 1:size(RemovedRegions,1)
        %SnRemRegionI = find(and(periods >= RemovedRegions(j,2) ,periods <= RemovedRegions(j,1) ));
         Sn(and(periods >= RemovedRegions(j,1) ,periods <= RemovedRegions(j,2) )) = 0;
    end
    
end

% plot(periods,Sn);
% jbase = find(periods >= Tp,1,'last');
[Energy,jj] = max(Sn);
% jj = jj + jbase -1;


TpC = periods(jj);
% plot([TpC,TpC],[0,1.2*Energy],'--k');

%Find turning points close to target
%First time it turns around
%To Right of TPCI in period space
a=Sn(1:jj-2);
b=Sn(2:jj-1);
c=Sn(3:jj);
TpCILB  = find((b<a & b<c)| (a ==0 | b ==0 | c == 0),1,'last')+1;

if isempty(TpCILB)
    TpCILB = 1;
end


%To Right of TPCI in period space
a=Sn(jj:end-2);
b=Sn(jj+1:end-1);
c=Sn(jj+2:end);
TpCIRB  = find((b<a & b<c)| (a ==0 | b ==0 | c == 0),1,'first')+jj;

if isempty(TpCIRB)
    TpCIRB = length(Sn) ;
end

clear a b c;

TpRegion = [periods(TpCIRB),periods(TpCILB)];
return


function [TpC,TpRegion,Energy] = fn_ExtractFirstPeak(Sn,periods)

% jbase = find(periods >= Tp,1,'last');
[Energy,jj] = max(Sn);
% jj = jj + jbase;

TpC = periods(jj);
TpCI = jj;

% plot([TpC,TpC],[0,1.2*Energy],'--r');

%Find turning points close to target
%First time it turns around
%To Right of TPCI in period space
a=Sn(1:TpCI-2);
b=Sn(2:TpCI-1);
c=Sn(3:TpCI);
TpCILB  = find(b<a & b<c ,1,'last')+1;

if isempty(TpCILB)
    TpCILB = 1;
end


%To Right of TPCI in period space
a=Sn(TpCI:end-2);
b=Sn(TpCI+1:end-1);
c=Sn(TpCI+2:end);
TpCIRB  = find(b<a & b<c ,1,'first')+TpCI;

if isempty(TpCIRB)
    TpCIRB = length(Sn) ;
end

clear a b c;


TpRegion = [periods(TpCIRB),periods(TpCILB)];


return





function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return


