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

function fn_ReadPlotAllSpectra(conc,TestName,Tp,probes)

if ~exist('conc','var') conc=1; end %39,79 or empty (1)
if ~exist('probes','var'); probes= 1:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var');     TestName= '17' ;end; %'14' '2','14'

MainFig = figure();

probeh = ceil(length(probes)/2);

%MainAxis  = subplot(probeh,2);
close all;
i = 1;
t0description ='waves reach x';
t1description ='final waves reach x';
% t1description = 'beach ref waves reach x';
for probei = probes
 %Read Data in Calib 1, at probe
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probei);
    
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
    [Sn,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
    period = 1./ff;
    
    %Get JonSwap
    SnJS  = jonswap(2*pi*ff,'wp',2*pi/Tp);
    SnJS = SnJS*max(Sn(period<2*Tp));
    
    SnS{i} = Sn(period<2*Tp);
    SnJSS{i} = SnJS(period<2*Tp);
    Times{i} = [round(tm(SI));round(tm(EI))];
    Periods{i} = period(period<2*Tp);
    i = i +1;
end

i = 1;
maxSn = max(cellfun(@max,SnS));
while i <= size(SnS,2)
    if i == 1
        
        if conc == 1
            Description = ['Calibration -  ',WaveType, ' Waves  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)', '   Between Times ',t0description,' [period = ',num2str(TLB),'] || ' ,t1description,' [period = ',num2str(TUB),']'  ];
        else
            Description = ['Conc = ',num2str(conc),' - ',WaveType, ' Waves  ', 'Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)', '   Between Times ' ,t0description,' [period = ',num2str(TLB),'] || ' ,t1description,' [period = ',num2str(TUB),']'  ];
        end
  
        sgtitle(Description);
        
    end
    subplot(2,probeh,i);
    plot(Periods{i},SnS{i},'-b','LineWidth',2);
    hold on;
    plot(Periods{i},SnJSS{i},'--k','LineWidth',2);
    xlabel('period (s)');
    ylabel('Spectra')
    xlim([0, 2*Tp]);
    ylim([0, 1.1*maxSn]);
    %xlim([0.6, 1.2]);
    title({['Probe ' , num2str(probes(i))], [ '(',num2str(Times{i}(1)),'[s],',num2str(Times{i}(2)),'[s])']});

    i = i+1;
end
return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

