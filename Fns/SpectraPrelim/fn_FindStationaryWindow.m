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

function [t0,t1] = fn_FindStationaryWindow(mean_times, WaveEnergyAtTargetPeriod,Tp,X,conc)

WaveMaker = 50;

% t0 = 50;
% t1 = mean_times(end);
% 
% tt=find(mean_times>=t0);

% WindowLength = 400;
% AVGmeantimes = zeros(1,length(mean_times)-2*WindowLength);
% AVGWETP = zeros(1,length(mean_times)-2*WindowLength);
% for i = WindowLength + 1:length(mean_times)-WindowLength
%     AVGmeantimes(i) = mean(mean_times(i - WindowLength : i + WindowLength));
%     AVGWETP(i) = mean(WaveEnergyAtTargetPeriod(i - WindowLength : i + WindowLength));
% end
% 
% WindowLengthN = 400;
% NewLength = floor(length(mean_times) / WindowLengthN);
% AVGmeantimes1 = zeros(1,NewLength);
% AVGWETP1 = zeros(1,NewLength);
% for i = 1:NewLength
%     AVGmeantimes1(i) = mean(mean_times((i-1)*WindowLengthN +1:(i)*WindowLengthN));
%     AVGWETP1(i) =  mean(WaveEnergyAtTargetPeriod((i-1)*WindowLengthN +1: (i)*WindowLengthN));
% end


if (conc == 39 || conc == 79)
    Tindic = fn_TestTimes(1.0/Tp,X,'attn');
else
    Tindic = fn_TestTimes(1.0/Tp,X,'calibration');
end

% if conc == 39 || conc == 79
%     t0 = Tindic(7).time;
%     t1 = Tindic(8).time;
% else
%     t0 = Tindic(7).time;
%     t1 = Tindic(9).time;   
% end
t0 = Tindic(7).time;
t1 = Tindic(8).time;

TimeBetweenArrivalAndReflection = t1 - t0;




figure();
plot(mean_times(:),WaveEnergyAtTargetPeriod(:),'-b');
hold on;
% plot([mean_times(1) , mean_times(end) ],[0,0],'--k');
% plot([mean_times(1) , mean_times(end) ],[10*WaveEnergyAtTargetPeriod(tt(1)),10*WaveEnergyAtTargetPeriod(tt(1))],'--k');
% plot(AVGmeantimes,AVGWETP,'-r');
% plot(AVGmeantimes1,AVGWETP1,'-g');

plot([t0 , t0 ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--k');
plot([t1 , t1 ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--k');


Tpers = fn_Tpers(Tp);
Tind = Tindic;
tvec=[]; count=1;
for loop=1:length(Tind)
if strfind(Tind(loop).description,'x')
tvec(count)=Tind(loop).time;
tvecdes{count}=Tind(loop).description; count=count+1;
end
end % end loop Tind
clear count Tind
[t0,it0]=min(tvec); tvec(it0)=[]; tvecdes(it0)=[]; clear it0
[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);
clear tvec

plot([t0 , t0 ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--r');
plot([t1 , t1 ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--r');

% plot([t0 , t0 ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--k');
% plot([t1 , t1 ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--k');
% title('Wave Energy vs Mean Time Of Window');
% 
% figure();
% plot(mean_times(tt(2:end)),WaveEnergyAtTargetPeriod(tt(2:end)) - WaveEnergyAtTargetPeriod(tt(1:end-1)),'-r');
% hold on;
% title('Wave Energy vs Mean Time Of Window');


%bisection method with new std measure
%might be interesting for something else,  but not so useful here
% CurrMean = mean(WaveEnergyAtTargetPeriod(:));
% plot([mean_times(1) , mean_times(end) ],[CurrMean,CurrMean],'--k');
% 
% CurrX0I = 1;
% CurrX1I = length(mean_times);
% for j = 1:3
%     
%     BisecI = floor((CurrX0I + CurrX1I)/2);
%     for i =CurrX0I: BisecI
%         lengthLHS =  BisecI - CurrX0I;
%         %sumsqLHS = sqrt((1.0/lengthLHS)*sum((WaveEnergyAtTargetPeriod(CurrX0I: BisecI) - CurrMean ).^2));
%         sumsqLHS = std(WaveEnergyAtTargetPeriod(CurrX0I: BisecI));
%     end
%     
%     for i = BisecI:CurrX1I
%         lengthRHS = CurrX1I -  BisecI;
%         %sumsqRHS = sqrt((1.0/lengthRHS)*sum((WaveEnergyAtTargetPeriod( BisecI:CurrX1I) - CurrMean ).^2));
%         sumsqRHS = std(WaveEnergyAtTargetPeriod(BisecI:CurrX1I));
%     end
%     
%     if sumsqLHS < sumsqRHS
%        CurrX0I = CurrX0I;
%        CurrX1I = BisecI;
%        CurrMean = mean(WaveEnergyAtTargetPeriod(CurrX0I: CurrX1I));
%     else
%        CurrX0I = BisecI;
%        CurrX1I = CurrX1I;  
%        CurrMean = mean(WaveEnergyAtTargetPeriod(CurrX0I: CurrX1I));
%        plot([mean_times(1) , mean_times(end) ],[CurrMean,CurrMean],'--k');
%     end
% end
% 
% plot([mean_times(CurrX0I) , mean_times(CurrX0I) ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--k');
% plot([mean_times(CurrX1I) , mean_times(CurrX1I) ],[0,1.1*max(WaveEnergyAtTargetPeriod)],'--k');
    


return



