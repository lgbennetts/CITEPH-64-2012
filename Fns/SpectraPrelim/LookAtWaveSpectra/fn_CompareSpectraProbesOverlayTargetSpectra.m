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

function fn_CompareSpectraCalibConcOverlayIndividual(conc,TestName,Tp,probes)

if ~exist('conc','var') conc=1; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= 11:20 ;end;     %probes= 1:20 ;end; 

if ~exist('probeiR','var'); probeiR= 1 ;end; 
if ~exist('probeiT','var'); probeiT= 2 ;end;%probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= '1'; end ;   %TestName= {'2','14'},{'8','16'} {'12','17'}


if ~exist('LBFactor','var'); LBFactor= 1 ;end; 
if ~exist('UBFactor','var'); UBFactor= 1 ;end; 
%{'14','9'}

%{'15','10'}
%{'16','8'}
% close all;

%MainAxis  = subplot(probeh,2);
i = 1;
t0description ='waves reach x';
t1description ='final waves reach x';

 %Read Data in Calib 1, at probe
[ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probeiR);

Tp = c_pram.period/10.0;
Hs = 10*c_pram.wave_height;

dt = tm(2) - tm(1);

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower

TpLB = Tp*LBFactor;
TpUB = Tp*UBFactor;

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower
TindLB = fn_Tind(conc,TpLB,ProbeLocXY(1),WaveType);
TindUB = fn_Tind(conc,TpUB,ProbeLocXY(1),WaveType);


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

% Periods1 = period(and(period> TpLB,period< TpUB));
% 
% SnS1 = Sn1(and(period> TpLB,period< TpUB));
% Times1 = [round(tm(SI));round(tm(EI))];

Periods1 = period;
SnS1 = sqrt(2*Sn1);
% SnS1 = Sn1;
Times1 = [round(tm(SI));round(tm(EI))];



[ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probeiT);

Tp = c_pram.period/10.0;
Hs = 10*c_pram.wave_height;

dt = tm(2) - tm(1);

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower

TpLB = Tp*LBFactor;
TpUB = Tp*UBFactor;

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower
TindLB = fn_Tind(conc,TpLB,ProbeLocXY(1),WaveType);
TindUB = fn_Tind(conc,TpUB,ProbeLocXY(1),WaveType);


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


% Periods2 = period(and(period> TpLB,period< TpUB));
% 
% SnS2 = Sn2(and(period> TpLB,period< TpUB));

%Convert to amplitude
Periods2 = period;
SnS2 = Sn2;
SnS2 = sqrt(2*Sn2);
Times2 = [round(tm(SI));round(tm(EI))];

JSspec  = jonswap(2*pi*ff,'wp',2*pi/Tp,'Hs',Hs/1000);
% JSspec = sqrt(2*JSspec);


figure();
hold on;
maxSn = max(max(SnS2),max(SnS1) );

t0str = num2str(Times1(1));
t1str = num2str(Times1(2));
Lab = [num2str(conc) '% Concentration ', WaveType ,' Probe ', num2str(probeiR),' (min:', t0str,' [s], ',t1str ' [s])'];
ipR1 = plot(Periods1,SnS1,'-r','DisplayName',Lab,'LineWidth',1);       

t0str = num2str(Times1(1));
t1str = num2str(Times1(2));
Lab = [num2str(conc) '% Concentration ', WaveType ,' Probe ', num2str(probeiT),' (min:', t0str,' [s], ',t1str ' [s])'];
ipR2 = plot(Periods2,SnS2,'-b','DisplayName',Lab,'LineWidth',1);       
   

t0str = num2str(Times1(1));
t1str = num2str(Times1(2));
 Description = {['Comparison of Individual Probes (',num2str(probeiR),' , ',num2str(probeiT),')', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)'],...
   ['Fourier Time Window - [' t0description,'( ',t0str,' [s])',' || ', t1description, '( ',t1str,' [s])']  };
title(Description);

%Jonswap
% plot(Periods2,JSspec,'--k','LineWidth',2);

% xlim([0, 2*Tp]);
% ylim([0, maxSn]);
xlabel('period (s)');
ylabel('Amplitude');
legend([ipR1,ipR2]);


% figure();
% hold on;
% Description = 'Basic Transmission Coefficients Plot';
% title(Description);
% plot(Periods1,movmean(InitialTimeAverage2,6) ./ movmean(InitialTimeAverage1,6),'.k');
% xlim([0, 2*Tp]);
% ylim([0, 1]);
% xlabel('period (s)');
% ylabel('Transmission Coefficients');
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

