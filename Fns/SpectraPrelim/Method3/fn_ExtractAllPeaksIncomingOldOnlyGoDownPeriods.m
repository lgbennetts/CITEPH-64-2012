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

function [Tps,TpsRegions] = fn_ExtractAllPeaksIncoming(conc,TestName,probes,Tp)


    [Sn,periods] = fn_ReadandCombineSpectra(conc,TestName,Tp,probes);
    
    figure();
    plot(periods,Sn,'-b');
    hold on;
    
    [TpC , TpRegion,PeakEnergy ] = fn_ExtractFirstPeak(Sn,periods);
    
    i= 1 ;
    Tps(i) = TpC;
    TpsRegions(i,:) = TpRegion ;
    Energy = PeakEnergy;
    while and(Energy > PeakEnergy/100.0, i < 20)
        i = i +1;
        [TpC , TpRegion,Energy ] = fn_ExtractOtherPeaks(Sn,periods,TpsRegions(i-1,2));
        Tps(i) = TpC;
        TpsRegions(i,:) = TpRegion ;
    end
    
    


return

function [TpC,TpRegion,Energy] = fn_ExtractOtherPeaks(Sn,periods,Tp)
jbase = find(periods >= Tp,1,'last');
[Energy,jj] = max(Sn(jbase:end));
jj = jj + jbase -1;


TpC = periods(jj);
plot([TpC,TpC],[0,1.2*Energy],'--k');

%Find turning points close to target
%First time it turns around
%To Right of TPCI in period space
a=Sn(1:jj-2);
b=Sn(2:jj-1);
c=Sn(3:jj);
TpCILB  = find(b<a & b<c,1,'last')+1;

if isempty(TpCILB)
    TpCILB = 1;
end


%To Right of TPCI in period space
a=Sn(jj:end-2);
b=Sn(jj+1:end-1);
c=Sn(jj+2:end);
TpCIRB  = find(b<a & b<c,1,'first')+jj;

if isempty(TpCIRB)
    TpCIRB = length(Sn) ;
end

clear a b c;

TpRegion = [periods(TpCILB),periods(TpCIRB)];
return


function [TpC,TpRegion,Energy] = fn_ExtractFirstPeak(Sn,periods)

% jbase = find(periods >= Tp,1,'last');
[Energy,jj] = max(Sn);
% jj = jj + jbase;

TpC = periods(jj);
TpCI = jj;

plot([TpC,TpC],[0,1.2*Energy],'--k');

%Find turning points close to target
%First time it turns around
%To Right of TPCI in period space
a=Sn(1:TpCI-2);
b=Sn(2:TpCI-1);
c=Sn(3:TpCI);
TpCILB  = find(b<a & b<c | (a ==0 | b ==0 | c == 0),1,'last')+1;

if isempty(TpCILB)
    TpCILB = 1;
end


%To Right of TPCI in period space
a=Sn(TpCI:end-2);
b=Sn(TpCI+1:end-1);
c=Sn(TpCI+2:end);
TpCIRB  = find(b<a & b<c | (a ==0 | b ==0 | c == 0),1,'first')+TpCI;

if isempty(TpCIRB)
    TpCIRB = length(Sn) ;
end

clear a b c;


TpRegion = [periods(TpCILB),periods(TpCIRB)];


return





function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return


