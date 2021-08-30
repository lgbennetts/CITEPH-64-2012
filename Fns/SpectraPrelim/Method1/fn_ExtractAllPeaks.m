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

function [Energy,Peaks] = fn_ExtractAllPeaks(tm,disp,Tp,conc,PX,dt)

    %first one
    [TgSpec,TgPer,nxtTp,DataProc] = fn_ExtractOnePeak(tm,disp,Tp,conc,PX,dt);
    Energy = [sum(TgSpec)];
    [~,TpCI] = max(TgSpec);
    Peaks = [TgPer(TpCI)];
    
    
    [TgSpec1,TgPer1,nxtTp1,DataProc1] = fn_ExtractOnePeak(tm,DataProc,nxtTp,conc,PX,dt);
    %Remove energy extracted in the previous one
    
    
%     [TgSpec2,TgPer2,nxtTp2,DataProc2] = fn_ExtractOnePeak(tm,DataProc1,nxtTp1,conc,PX,dt);

%


return


function [TgSpec,TgPer,nxtTp,DataProc] = fn_ExtractOnePeak(tm,disp,Tp,conc,PX,dt)
% time between waves reach and MIZ reflection
%Bound the speed, lets pick a window
Tind = fn_Tind(conc,Tp,PX);
TindA = fn_Tind(conc,Tp*0.99,PX); %Slower
TindB = fn_Tind(conc,Tp*1.01,PX); % Faster

%pick t0 as slowest possible time
%pick t1 as fastest posssible time

t0 = TindA(7).time;
t1 = TindB(8).time;


dtlength = (t1 - t0) / dt;

timewindowlength = 2^floor(log2(dtlength));
SI = floor((t0 -tm(1)) / dt);
EI = SI + timewindowlength-1;

%% Fourier Transform
% TwindowIndex = find(and(tm <= t1, tm >= t0));
% 
[Sn,an,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
periods = 1./ff;

% [Peak,TpCI] = max(Sn);
% TpC = periods(TpCI);

%Find location close to target period
[~,TpCI]=min(abs(periods-Tp));
Peak = Sn(TpCI);
TpC = periods(TpCI);


%Find turning points close to target
%First time it turns around
%To Right of TPCI in period space
a=Sn(1:TpCI-2);
b=Sn(2:TpCI-1);
c=Sn(3:TpCI);
TpCILB  = find(b<a & b<c,1,'last')+1;

if isempty(TpCILB)
    TpCILB = 1;
end


%To Right of TPCI in period space
a=Sn(TpCI:end-2);
b=Sn(TpCI+1:end-1);
c=Sn(TpCI+2:end);
TpCIRB  = find(b<a & b<c,1,'first')+TpCI;

if isempty(TpCIRB)
    TpCIRB = length(Sn) ;
end

clear a b c;

%find closest peak in the region around target period
[Peak,TpCI] = max(Sn(TpCILB:TpCIRB));
TpCI = TpCILB + TpCI -1;
TpC = periods(TpCI);

%Find turning points close to target
%First time it turns around
%To Right of TPCI in period space
a=Sn(1:TpCI-2);
b=Sn(2:TpCI-1);
c=Sn(3:TpCI);
TpCILB  = find(b<a & b<c,1,'last')+1;

if isempty(TpCILB)
    TpCILB = 1;
end


%To Right of TPCI in period space
a=Sn(TpCI:end-2);
b=Sn(TpCI+1:end-1);
c=Sn(TpCI+2:end);
TpCIRB  = find(b<a & b<c,1,'first')+TpCI;

if isempty(TpCIRB)
    TpCIRB = length(Sn) ;
end

clear a b c;



figure();
plot(periods,Sn,'-');
hold on;
plot([periods(TpCIRB),periods(TpCIRB)],[0,Peak],'--k');
plot([periods(TpCILB),periods(TpCILB)],[0,Peak],'--k');

TgSpec = Sn(TpCILB :TpCIRB);
TgPer = periods(TpCILB :TpCIRB);

Sn(TpCILB :TpCIRB) = zeros(1,TpCIRB-TpCILB+1);

[~,nxtTpI] = max(Sn);
nxtTp = periods(nxtTpI);

halffreq = length(an)/2;

ancurr = zeros(size(an));

ancurr(end -TpCIRB : end -TpCILB ) = an(end -TpCIRB : end -TpCILB );
ancurr(TpCILB : TpCIRB ) = an(TpCILB : TpCIRB );

figure();
plot(tm,disp,'-b');
hold on;

meandisp = mean(disp((SI:EI)));
dispA =fft(ancurr);
perdispW = meandisp+ real(dispA);

ShiftToZeroT = rem(SI,length(perdispW));
perdispWAt0 = circshift(perdispW,ShiftToZeroT);
CopyArrayDisp = repmat(perdispWAt0, ceil(length(tm)/length(perdispWAt0)));
CopyArrayDisp = CopyArrayDisp(1 : length(tm));
CopyArrayDisp = CopyArrayDisp.';


% perdispAll = ;

%extend signal
% 
% plot(tm,CopyArrayDisp,'-r');
% 
% 
plot(tm(SI:EI),disp(SI:EI),'--k');
plot(tm(SI:EI),perdispW,'--g');
% 
% 
% plot(tm,disp - CopyArrayDisp,'--k');

DataProc = disp - CopyArrayDisp;

plot(tm,DataProc,'-r');
return





function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return


