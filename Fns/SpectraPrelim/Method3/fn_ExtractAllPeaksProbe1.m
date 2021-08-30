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

function [Energy,Tps,Regions,TgTps] = fn_ExtractAllPeaksProbe1(tm,disp,Tp,conc,PX,dt)

    %first one
    RemovedRegions = [100*Tp,min(4,3*Tp); Tp/10.0,0];
    TgTps = [Tp];
    [EnergyAtTp, RemoveRegion,NxtTp,TpC] = fn_ExtractOnePeak(tm,disp,Tp,conc,PX,dt,RemovedRegions);
    
    Energy = [EnergyAtTp];
    Tps = [TpC];
    RemovedRegions = [RemovedRegions;RemoveRegion];
    
    while and(Energy(end) >   Energy(1)/100, length(Energy) < 5 )
        TgTps = [TgTps;NxtTp];
        [EnergyAtTp, RemoveRegion,NxtTp,TpC] =  fn_ExtractOnePeak(tm,disp,NxtTp,conc,PX,dt,RemovedRegions);
        Energy = [Energy;EnergyAtTp];
        Tps = [Tps;TpC];
        RemovedRegions = [RemovedRegions;RemoveRegion];
    end
    
    Regions = RemovedRegions(3:end,:);

    
%     Energy = [sum(TgSpec)];
%     [~,TpCI] = max(TgSpec);
%     Peaks = [TgPer(TpCI)];
%     
%     
%     [TgSpec1,TgPer1,nxtTp1,DataProc1] = fn_ExtractOnePeak(tm,DataProc,,conc,PX,dt);
    %Remove energy extracted in the previous one
    
    
%     [TgSpec2,TgPer2,nxtTp2,DataProc2] = fn_ExtractOnePeak(tm,DataProc1,nxtTp1,conc,PX,dt);

%


return


function [EnergyAtTp, RemoveRegion,NxtTp,TpC] = fn_ExtractOnePeak(tm,disp,Tp,conc,PX,dt,RemovedRegions)
% time between waves reach and MIZ reflection
%Bound the speed, lets pick a window
Tind = fn_Tind(conc,Tp,PX);
TindA = fn_Tind(conc,Tp*0.95,PX); %Slower
TindB = fn_Tind(conc,Tp*1.05,PX); % Faster

%pick t0 as slowest possible time
%pick t1 as fastest posssible time

t0 = TindA(7).time;
t1 = TindB(8).time;
% t0 = Tind(7).time;
% t1 = Tind(8).time;


dtlength = (t1 - t0) / dt;

timewindowlength = 2^floor(log2(dtlength));
SI = floor((t0 -tm(1)) / dt);
EI = SI + timewindowlength-1;

%% Fourier Transform
[Sn,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
periods = 1./ff;

%Zero out remove regions
if ~isempty(RemovedRegions)
    for j = 1:size(RemovedRegions,1)
%         SnRemRegionI = find(and(periods >= RemovedRegions(j,2) ,periods <= RemovedRegions(j,1) ));
        SNIRB = find(periods > RemovedRegions(j,2),1,'last');
        SNILB = find(periods < RemovedRegions(j,1),1,'first');
        Sn(SNILB : SNIRB) = zeros(1,SNIRB -SNILB +1);
%         Sn(min(SnRemRegionI +1 , length(Sn) )) = 0;
    end
    
end

%Find peak 
[Peak,TpCI]=max(Sn);
TpC = periods(TpCI);


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


RemoveRegion = [periods(TpCILB),periods(TpCIRB)];

% figure();
% plot(periods,Sn,'-');
% hold on;
% plot([periods(TpCIRB),periods(TpCIRB)],[0,Peak],'--k');
% plot([periods(TpCILB),periods(TpCILB)],[0,Peak],'--k');
% hold off;

% TgSpec = Sn(TpCILB :TpCIRB);
EnergyAtTp = sum(Sn(TpCILB :TpCIRB));
% TgPer = periods(TpCILB :TpCIRB);

Sn(TpCILB :TpCIRB) = zeros(1,TpCIRB-TpCILB+1);

[~,NxtTpI] = max(Sn);
NxtTp = periods(NxtTpI);


return





function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return


