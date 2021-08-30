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

function fn_ReadPlotAllRegIrregSpectra(conc,TestName,Tp,probes)

if ~exist('conc','var') conc=1; end %39,79 or empty (1)
if ~exist('probes','var'); probes= 11:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= '14'; end ;   %TestName= {'2','14'},{'8','16'} {'12','17'}

close all;
MainFig = figure();
hold on;

probeh = ceil(length(probes)/2);

%MainAxis  = subplot(probeh,2);
i = 1;
t0descriptionr ='waves reach x';
t0descriptionb ='beach ref waves reach x';

t1descriptionf ='final waves reach x';
t1descriptionb ='beach ref waves reach x';
% t1description ='reflected waves from beach reach wave maker';
% t1description = 'beach ref waves reach x';
for probei = probes
 %Read Data in Calib 1, at probe
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probei);
    
    Tp = c_pram.period/10.0;
    Hs = 10*c_pram.wave_height;
    
    dt = tm(2) - tm(1);

    Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower
    
    TLB = Tp;
    TUB = Tp;
    
    TindLB = fn_Tind(conc,TLB,ProbeLocXY(1),WaveType);
    TindUB = fn_Tind(conc,TUB,ProbeLocXY(1),WaveType);

    t0index = find(strcmp({Tind.description},t0descriptionr));
    t0 = TindLB(t0index).time;
    t1index = find(strcmp({Tind.description},t1descriptionb));
    t1 = TindUB(t1index).time ;
    
    dtlength = (t1 - t0) / dt;
    timewindowlength = 2^floor(log2(dtlength));
    
    
    MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
    SI = MI - timewindowlength/2;
    EI = MI + timewindowlength/2 -1;

    %% Fourier Transform
    [Sn1,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;

    Periods1{i} = period(period<2*Tp);
    
    SnS1{i} = Sn1(period<2*Tp);
    Times1{i} = [round(tm(SI));round(tm(EI))];
    
    
    %second on
    timebetween = t1 - t0;
    t0index = find(strcmp({Tind.description},t0descriptionb));
    t0 = TindLB(t0index).time; 
    t1 = t0 + timebetween;
    
    dtlength = (t1 - t0) / dt;
    timewindowlength = 2^floor(log2(dtlength));

    MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
    SI = MI - timewindowlength/2;
    EI = MI + timewindowlength/2 -1;

    
    %% Fourier Transform
    [Sn2,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;

    Periods2{i} = period(period<2*Tp);
    
    SnS2{i} = Sn2(period<2*Tp);
    Times2{i} = [round(tm(SI));round(tm(EI))];
    
    i = i +1;
end

i = 1;
maxSn = max(max(cellfun(@max,SnS2)),max(cellfun(@max,SnS1)) );
while i <= size(SnS1,2)

        if i == 1
            t0str = num2str(max(cellfun(@min,Times1(1:10))));
            t1str = num2str(min(cellfun(@max,Times1(1:10))));
            Lab = ['Irregular Reflection Probes (min:', t0str,' [s], ',t1str ' [s])'];
            ipR1 = plot(Periods1{i},SnS1{i},'-r','DisplayName',Lab,'LineWidth',1,'Color','#EDB120');       
        else
            plot(Periods1{i},SnS1{i},'-r','LineWidth',1,'Color','#EDB120'); 
        end
            
        if i == 1
            t0str = num2str(max(cellfun(@min,Times2(1:10))));
            t1str = num2str(min(cellfun(@max,Times2(1:10))));
            Lab = ['Irregular Reflection Probes (min:', t0str,' [s], ',t1str ' [s])'];
            ipR2 = plot(Periods2{i},SnS2{i},'-b','DisplayName',Lab,'LineWidth',1,'Color',	'#4DBEEE');       
        else
            plot(Periods2{i},SnS2{i},'-b','LineWidth',1,'Color',	'#4DBEEE'); 
        end     
    i = i+1;
end

SnS1L = [SnS1{1},SnS1{2},SnS1{3},SnS1{4},SnS1{5},SnS1{6},SnS1{7},SnS1{8},SnS1{9},SnS1{10}];
SnS2L = [SnS2{1},SnS2{2},SnS2{3},SnS2{4},SnS2{5},SnS2{6},SnS1{7},SnS2{8},SnS2{9},SnS2{10}];
InitialTimeAverage1 = mean(SnS1L,2);
InitialTimeAverage2 = mean(SnS2L,2);



 Description = ['Calibration - Irregular Waves Compare Reflection Probes Different Time Windows  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)'];
title(Description);

Lab = ['Irregular Reflection Probes', t0descriptionr ' || ', t1descriptionb];
aR1 = plot(Periods1{1},InitialTimeAverage1,'--r','LineWidth',2,'DisplayName',Lab);

Lab = ['Irregular Reflection Probes', t0descriptionb ' ||  Same Time Length as Above'];
aR2 = plot(Periods2{1},InitialTimeAverage2,'--b','LineWidth',2,'DisplayName',Lab);

xlim([0, 2*Tp]);
ylim([0, maxSn]);
xlabel('period (s)');
ylabel('Spectra');
legend([ipR1,ipR2,aR1,aR2]);
hold off;

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

