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

function [Energys] = fn_EnergyAtPeaks(tm,disp,conc,PX,dt,Peaks_Baseline,Regions_Baseline,TgTps_Baseline)

       [Energy] = fn_ExtractOnePeak(tm,disp,TgTps_Baseline(1),Regions_Baseline(1,:),[],conc,PX,dt);
       Energys = [Energy];
       for i = 2: size(Peaks_Baseline,1)
           [Energy] = fn_ExtractOnePeak(tm,disp,TgTps_Baseline(i),Regions_Baseline(i,:),Regions_Baseline(1:i-1,:),conc,PX,dt);
           Energys = [Energys;Energy];
       end

return


function [EnergyAtTp] = fn_ExtractOnePeak(tm,disp,Tp,Regions,ZeroRegions,conc,PX,dt)
% time between waves reach and MIZ reflection
%Bound the speed, lets pick a window
Tind = fn_Tind(conc,Tp,PX);
TindA = fn_Tind(conc,Tp*0.95,PX); %Slower
TindB = fn_Tind(conc,Tp*1.05,PX); % Faster

%pick t0 as slowest possible time
%pick t1 as fastest posssible time

t0 = TindA(7).time;
t1 = TindB(8).time;


dtlength = (t1 - t0) / dt;

timewindowlength = 2^floor(log2(dtlength));
SI = floor((t0 -tm(1)) / dt);
EI = SI + timewindowlength-1;

%% Fourier Transform
[Sn,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
periods = 1./ff;

%         SNIRB = find(periods > RemovedRegions(j,2),1,'last');
%         SNILB = find(periods < RemovedRegions(j,1),1,'first');

if ~isempty(ZeroRegions)
    for j = 1:size(ZeroRegions,1)
%         SnRemRegionI = find(and(periods >= RemovedRegions(j,2) ,periods <= RemovedRegions(j,1) ));
        SNIRB = find(periods >= ZeroRegions(j,2),1,'last');
        SNILB = find(periods <= ZeroRegions(j,1),1,'first');
        Sn(SNILB : SNIRB) = zeros(1,SNIRB -SNILB +1);
%         Sn(min(SnRemRegionI +1 , length(Sn) )) = 0;
    end
    
end

TpCILB =  find(periods == Regions(1)) ;
TpCIRB =  find(periods == Regions(2)) ;

if isempty(TpCILB)
    TpCILB =  find(periods < Regions(1),1,'first') ;
end

if isempty(TpCIRB)
    TpCIRB = find(periods > Regions(2),1,'last');
end



if TpCILB == TpCIRB
    TpCILB = max(1,TpCILB -1);
    TpCIRB = min(length(Sn),TpCIRB + 1);
end

figure();
plot(periods,Sn,'-');
hold on;

plot([periods(TpCILB),periods(TpCILB)],[0,max(Sn)],'--k');
plot([periods(TpCIRB),periods(TpCIRB)],[0,max(Sn)],'--k');

hold off;

EnergyAtTp = sum(Sn(TpCILB :TpCIRB));


return





function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return


