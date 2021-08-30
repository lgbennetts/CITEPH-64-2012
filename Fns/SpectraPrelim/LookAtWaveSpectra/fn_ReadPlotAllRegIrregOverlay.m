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
if ~exist('probes','var'); probes= 1:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {'2','14'}; end ;   %TestName= {'2','14'},{'8','16'} {'12','17'}

close all;
MainFig = figure();
hold on;

probeh = ceil(length(probes)/2);

%MainAxis  = subplot(probeh,2);
i = 1;
t0description ='waves reach x';
t1description ='beach ref waves reach x';
% t1description ='reflected waves from beach reach wave maker';
% t1description = 'beach ref waves reach x';
for probei = probes
 %Read Data in Calib 1, at probe
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName{1}, probei);
    
    Tp = c_pram.period/10.0;
    Hs = 10*c_pram.wave_height;
    
    dt = tm(2) - tm(1);

    Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower
    %TindB = fn_Tind(conc,Tp,ProbeLocXY(1));
    
    TLB = Tp;
    TUB = Tp;
    
    TindLB = fn_Tind(conc,TLB,ProbeLocXY(1),WaveType);
    TindUB = fn_Tind(conc,TUB,ProbeLocXY(1),WaveType);
    
    t0index = find(strcmp({Tind.description},t0description));
    t0 = TindLB(t0index).time;
    t1index = find(strcmp({Tind.description},t1description));
    t1 = TindUB(t1index).time ;
    
    
    dtlength = (t1 - t0) / dt;

    timewindowlength = 2^floor(log2(dtlength));
%     [num2str(probei),'   ',num2str(timewindowlength),'  ',num2str(t0),'s  to ',num2str(t1),'s']

      %fixed start
%     SI = floor((t0 -tm(1)) / dt)+1;
%     EI = SI + timewindowlength-1;
%     
      %fixed end
%     EI = floor((t1 -tm(1)) / dt);
%     SI = EI - timewindowlength-1;

    MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
    SI = MI - timewindowlength/2;
    EI = MI + timewindowlength/2 -1;

    % [timewindowlength,tm(SI),tm(EI)]

    %% Fourier Transform
    [Sn1,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;

    Periods1{i} = period(period<2*Tp);
    
    Sn1S{i} = Sn1(period<2*Tp);
    Times1{i} = [round(tm(SI));round(tm(EI))];
    
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName{2}, probei);
    
    Tp = c_pram.period/10.0;
    Hs = 10*c_pram.wave_height;
    
    dt = tm(2) - tm(1);

    Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower
    %TindB = fn_Tind(conc,Tp,ProbeLocXY(1));
    
    TLB = Tp;
    TUB = Tp;
    
    TindLB = fn_Tind(conc,TLB,ProbeLocXY(1),WaveType);
    TindUB = fn_Tind(conc,TUB,ProbeLocXY(1),WaveType);
    
    t0index = find(strcmp({Tind.description},t0description));
    t0 = TindLB(t0index).time;
    t1index = find(strcmp({Tind.description},t1description));
    t1 = TindUB(t1index).time;  
    
    
    dtlength = (t1 - t0) / dt;

    timewindowlength = 2^floor(log2(dtlength));
%     [num2str(probei),'   ',num2str(timewindowlength),'  ',num2str(t0),'s  to ',num2str(t1),'s']

      %fixed start
%     SI = floor((t0 -tm(1)) / dt)+1;
%     EI = SI + timewindowlength-1;
%     
      %fixed end
%     EI = floor((t1 -tm(1)) / dt);
%     SI = EI - timewindowlength-1;

    MI = floor((((t0 + t1) / 2) - tm(1))/dt  );
    SI = MI - timewindowlength/2;
    EI = MI + timewindowlength/2 -1;

    % [timewindowlength,tm(SI),tm(EI)]

    %% Fourier Transform
    [Sn2,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;
    
    %Get JonSwap
    SnJS  = jonswap(2*pi*ff,'wp',2*pi/Tp);
    SnJS = SnJS*max(Sn2(period<2*Tp));
    
    Sn2S{i} = Sn2(period<2*Tp);
    SnJSS{i} = SnJS(period<2*Tp);
    Times2{i} = [round(tm(SI));round(tm(EI))];
    Periods2{i} = period(period<2*Tp);
    i = i +1;
end

i = 1;
maxSn = max(cellfun(@max,Sn2S));
while i <= size(Sn2S,2)
%     if i == 1
%         
%         if conc == 1
%             Description = ['Calibration -  ', WaveType ' Waves  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)', '   Between Times ',t0description,' [period = ',num2str(TLB),'] || ' ,t1description,' [period = ',num2str(TUB),']'  ];
%         else
%             Description = ['Conc = ',num2str(conc),' -' WaveType ' Waves  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)', '   Between Times ' ,t0description,' [period = ',num2str(TLB),'] || ' ,t1description,' [period = ',num2str(TUB),']'  ];
%         end
%   
%         sgtitle(Description);
%         
%     end
    
    %seperate plots
%     subplot(2,probeh,i);
%     plot(Periods2{i},Sn2S{i},'-b','LineWidth',2);
%     hold on;
%     plot(Periods1{i},Sn1S{i},'-r','LineWidth',1);
%     plot(Periods2{i},SnJSS{i},'--k','LineWidth',2);
%     xlabel('period (s)');
%     ylabel('Spectra')
%     xlim([0, 2*Tp]);
%     ylim([0, 1.1*maxSn]);
%     %xlim([0.6, 1.2]);
%     title({['Probe ' , num2str(probes(i))], [ '(',num2str(Times{i}(1)),'[s],',num2str(Times{i}(2)),'[s])']});


    if probes(i) < 11

        if (probes(i) ==1)
            t0str = num2str(max(cellfun(@min,Times2(1:10))));
            t1str = num2str(min(cellfun(@max,Times2(1:10))));
            Lab = ['Irregular Reflection Probes (min:', t0str,' [s], ',t1str ' [s])'];
            ipR = plot(Periods2{i},Sn2S{i},'-','DisplayName',Lab,'LineWidth',1,'Color',	'#4DBEEE');

            t0str = num2str(max(cellfun(@min,Times1(1:10))));
            t1str = num2str(min(cellfun(@max,Times1(1:10))));
            Lab = ['Regular Reflection Probes (min:', t0str,' [s], ',t1str ' [s])'];
            rpR = plot(Periods1{i},Sn1S{i},'-','DisplayName',Lab,'LineWidth',1,'Color','#EDB120');       
        end
            plot(Periods2{i},Sn2S{i},'-','LineWidth',1,'Color',	'#4DBEEE');
            plot(Periods1{i},Sn1S{i},'-','LineWidth',1,'Color','#EDB120'); 
            
    else
        
        if (probes(i) ==11)
            t0str = num2str(max(cellfun(@min,Times2(11:20))));
            t1str = num2str(min(cellfun(@max,Times2(11:20))));
            Lab = ['Irregular Transmission Probes  (min:', t0str,' [s], ',t1str ' [s])'];
            ipT = plot(Periods2{i},Sn2S{i},'-b','DisplayName',Lab,'LineWidth',1);

            t0str = num2str(max(cellfun(@min,Times2(11:20))));
            t1str = num2str(min(cellfun(@max,Times2(11:20))));
            Lab = ['Regular Transmission Probes  (min:', t0str,' [s], ',t1str ' [s])'];
            rpT = plot(Periods1{i},Sn1S{i},'-r','DisplayName',Lab,'LineWidth',1);       
        end
            plot(Periods2{i},Sn2S{i},'-b','LineWidth',1);
            plot(Periods1{i},Sn1S{i},'-r','LineWidth',1);         

    
    end 
    i = i+1;
end

%Plot average of spectra as well
[SnA1,periodsA1,SnA1Times] = fn_ReadandCombineSpectra(conc,TestName{1},Tp,probes,t0description,t1description);
[SnA2,periodsA2,SnA2Times] = fn_ReadandCombineSpectra(conc,TestName{2},Tp,probes,t0description,t1description);


if conc == 1
    Description = ['Calibration -  ', ' Regular and Irregular Waves  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)', '   Between Times ',t0description,' [period = ',num2str(TLB),'] || ' ,t1description,' [period = ',num2str(TUB),']'  ];
else
    Description = ['Conc = ',num2str(conc),' - Regular and Irregular Waves  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)', '   Between Times ' ,t0description,' [period = ',num2str(TLB),'] || ' ,t1description,' [period = ',num2str(TUB),']'  ];
end
title(Description);

% figure();
% 
% plot(periodsA2,SnA2,'-b','LineWidth',2);
% hold on;
% plot(periodsA1,SnA1,'-r','LineWidth',1);
% 
% SnJSA2  = jonswap(2*pi*ff,'wp',2*pi/Tp);
% SnJSA2 = SnJSA2*max(SnA2(period<2*Tp));
% xlim([0, 2*Tp]);
% ylim([0, 1.1*max(SnA2)]);
% plot(periodsA2,SnJSA2,'--k','LineWidth',2);
% title(Description);
% xlabel('period (s)');
% ylabel('Spectra');
% hold off;

Lab = 'Average Of All Probes For Irregular Spectra';
as = plot(periodsA2,SnA2,'--k','DisplayName',Lab,'LineWidth',2);

Lab = 'JonSwap Spectrum';
SnJSA2  = jonswap(2*pi./periodsA2,'wp',2*pi/Tp);
SnJSA2 = SnJSA2*max(SnA2(periodsA2<2*Tp));

%SnJSA2 = SnJSA2*maxSn;
js = plot(periodsA2,SnJSA2,'--g','DisplayName',Lab,'LineWidth',2);

xlim([0, 2*Tp]);
ylim([0, maxSn]);
xlabel('period (s)');
ylabel('Spectra');
legend([ipR,rpR,ipT,rpT,as,js]);
hold off;


return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

