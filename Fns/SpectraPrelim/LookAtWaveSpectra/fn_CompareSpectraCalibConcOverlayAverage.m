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

function fn_CompareSpectraCalibConcOverlayAverage(conc,TestName,Tp,probes)

if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= 11:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {'14','9'}; end ;   %TestName= {'2','14'},{'8','16'} {'12','17'}
%{'14','9'}
%{'15','10'}
%{'16','8'}
% close all;

probeh = ceil(length(probes)/2);

%MainAxis  = subplot(probeh,2);
i = 1;
t0description ='waves reach x';
t1description ='final waves reach x';

for probei = probes
 %Read Data in Calib 1, at probe
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{1}, probei);
    
    Tp = c_pram.period/10.0;
    Hs = 10*c_pram.wave_height;
    
    dt = tm(2) - tm(1);

    TpLB = Tp*0.5;
    TpUB = Tp*2;
    
    Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
    
    TindLB = fn_Tind(conc{1},TpLB,ProbeLocXY(1),WaveType);
    TindUB = fn_Tind(conc{1},TpUB,ProbeLocXY(1),WaveType);

    t0index = find(strcmp({TindLB.description},t0description));
    t0 = TindLB(t0index).time;
    t1index = find(strcmp({TindUB.description},t1description));
    t1 = TindUB(t1index).time ;
    
    dtlength = (t1 - t0) / dt;
    timewindowlength = 2^floor(log2(dtlength));
    
    
    MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
    SI = MI - timewindowlength/2;
    EI = MI + timewindowlength/2 -1;

    %% Fourier Transform
    [Sn1,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;

    Periods1{i} = period(and(period>TpLB,period<TpUB));
    
    SnS1{i} = Sn1(and(period>TpLB,period<TpUB));
    Times1{i} = [round(tm(SI));round(tm(EI))];
    

    
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{2}, probei);
    
    Tp = c_pram.period/10.0;
    Hs = 10*c_pram.wave_height;
    
    dt = tm(2) - tm(1);

    TpLB = Tp*0.5;
    TpUB = Tp*2;
    
    Tind = fn_Tind(conc{2},Tp,ProbeLocXY(1),WaveType); %Slower
    
    TindLB = fn_Tind(conc{2},TpLB,ProbeLocXY(1),WaveType);
    TindUB = fn_Tind(conc{2},TpUB,ProbeLocXY(1),WaveType);

    t0index = find(strcmp({TindLB.description},t0description));
    t0 = TindLB(t0index).time;
    t1index = find(strcmp({TindUB.description},t1description));
    t1 = TindUB(t1index).time ;
    
    dtlength = (t1 - t0) / dt;
    timewindowlength = 2^floor(log2(dtlength));
    
    
    MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
    SI = MI - timewindowlength/2;
    EI = MI + timewindowlength/2 -1;

    %% Fourier Transform
    [Sn2,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;


    Periods2{i} = period(and(period>TpLB,period<TpUB));
    
    SnS2{i} = Sn2(and(period>TpLB,period<TpUB));
    Times2{i} = [round(tm(SI));round(tm(EI))];
    
    i = i +1;
end

i = 1;
figure();
hold on;
maxSn = max(max(cellfun(@max,SnS2)),max(cellfun(@max,SnS1)) );
while i <= size(SnS1,2)

        if i == 1
            t0str = num2str(max(cellfun(@min,Times1(1:10))));
            t1str = num2str(min(cellfun(@max,Times1(1:10))));
            Lab = ['No Disks Irregular Transmission Probes (min:', t0str,' [s], ',t1str ' [s])'];
            ipR1 = plot(Periods1{i},SnS1{i},'-r','DisplayName',Lab,'LineWidth',1);       
        else
            plot(Periods1{i},SnS1{i},'-r','LineWidth',1); 
        end
            
        if i == 1
            t0str = num2str(max(cellfun(@min,Times2(1:10))));
            t1str = num2str(min(cellfun(@max,Times2(1:10))));
            Lab = ['39% Concentration Irregular Transmission Probes (min:', t0str,' [s], ',t1str ' [s])'];
            ipR2 = plot(Periods2{i},SnS2{i},'-b','DisplayName',Lab,'LineWidth',1);       
        else
            plot(Periods2{i},SnS2{i},'-b','LineWidth',1); 
        end     
    i = i+1;
end

t0str = num2str(max(cellfun(@min,Times1(1:10))));
t1str = num2str(min(cellfun(@max,Times1(1:10))));
 Description = {['Individual Compare Disks to Calibation - Transmission Probes', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)'],...
   ['Fourier Time Window - [' t0description,'( ',t0str,' [s])',' || ', t1description, '( ',t1str,' [s])']  };
title(Description);

xlim([0, 2*Tp]);
ylim([0, maxSn]);
xlabel('period (s)');
ylabel('Spectra');
legend([ipR1,ipR2]);


SnS1L = [SnS1{1},SnS1{2},SnS1{3},SnS1{4},SnS1{5},SnS1{6},SnS1{7},SnS1{8},SnS1{9},SnS1{10}];
SnS2L = [SnS2{1},SnS2{2},SnS2{3},SnS2{4},SnS2{5},SnS2{6},SnS1{7},SnS2{8},SnS2{9},SnS2{10}];
InitialTimeAverage1 = mean(SnS1L,2);
InitialTimeAverage2 = mean(SnS2L,2);


figure();

Lab = ['No Disks Irregular Transmission Probes' ];
aR1 = plot(Periods1{1},InitialTimeAverage1,'-r','LineWidth',2,'DisplayName',Lab);
hold on;
Lab = ['39% Concentration Irregular Transmission Probes'];
aR2 = plot(Periods2{1},InitialTimeAverage2,'-b','LineWidth',2,'DisplayName',Lab);

Description = {['Average Compare Disks to Calibation - Transmission Probes', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)'],...
   ['Fourier Time Window - [' t0description,'( ',t0str,' [s])',' || ', t1description, '( ',t1str,' [s])']  };
title(Description);

xlim([0, 2*Tp]);
ylim([0, maxSn]);
xlabel('period (s)');
ylabel('Spectra');
legend([aR1,aR2]);
hold off;

% figure();
% hold on;
% Description = 'Basic Transmission Coefficients Plot';
% title(Description);
% plot(Periods1{1},movmean(InitialTimeAverage2,6) ./ movmean(InitialTimeAverage1,6),'.k');
% xlim([0, 2*Tp]);
% ylim([0, 1]);
% xlabel('period (s)');
% ylabel('Transmission Coefficients');


[N,edges,bin] = histcounts(Periods1{1},20);
ITA1 = accumarray(bin(:),InitialTimeAverage1,[],@sum);
ITA2 = accumarray(bin(:),InitialTimeAverage2,[],@sum);

EdgesH = (edges(2:end) + edges(1:end-1)) /2;

figure();
hold on;
Description = 'Binned Transmission Coefficients Plot';
title(Description);
plot(EdgesH,ITA2 ./ ITA1,'.k');
xlim([0, 2*Tp]);
ylim([0, 1]);
xlabel('period (s)');
ylabel('Transmission Coefficients');

% hold on;
% hold off;

% figure();
% hold on;
% title(Description);
% Lab = [' Average of Reflection Probes Between Times ',t0description,' [period = ',num2str(TLB),'] || ' ,t1descriptionb,' [period = ',num2str(TUB),']'  ];
% aR = plot(periodsR,SnR,'-r','DisplayName',Lab,'LineWidth',2);
% 
% Lab = [' Average of Transmission Probes Between Times ',t0description,' [period = ',num2str(TLB),'] || ' ,t1descriptionf,' [period = ',num2str(TUB),']'  ];
% aT = plot(periodsT,SnT,'-b','DisplayName',Lab,'LineWidth',2);
% 
% 
% Lab = [' Average of Transmission (Summed Into Lower Periods) Probes Between Times ',t0description,' [period = ',num2str(TLB),'] || ' ,t1descriptionf,' [period = ',num2str(TUB),']'  ];
% aTAR = plot(periodsR,SnTAR,'--k','DisplayName',Lab,'LineWidth',2);

% 
% Lab = 'JonSwap Spectrum';
% SnJSA2  = jonswap(2*pi./periodsA2,'wp',2*pi/Tp);
% SnJSA2 = SnJSA2*max(SnA2(periodsA2<2*Tp));
% 
% %SnJSA2 = SnJSA2*maxSn;
% js = plot(periodsA2,SnJSA2,'--g','DisplayName',Lab,'LineWidth',2);

% xlim([0, 2*Tp]);
% ylim([0, maxSn]);
% xlabel('period (s)');
% ylabel('Spectra');
% legend([ipR1,ipR2]);
% hold off;


return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

