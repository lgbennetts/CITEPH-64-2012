close all;
clear all;

% probesRef = 1:10;
% probesTra = 11:20;
probesRef = [1];
probesTra = [11];

CalibRef = [];
CalibRefTDiff = [];
ExpTra = [];
ExpTraTDiff = [];


%Tp = 0.8, Hs = 0.02
% TNameCalib = '14';
% TNameExp = '9';

% %Tp = 1.4, Hs = 0.04
% TNameCalib = '15';
% TNameExp = '10';

% %Tp = 1.4, Hs = 0.08
TNameCalib = '16';
TNameExp = '8';


%Calibration

for i = 1:size(probesTra,2)
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(1,TNameCalib,probesRef(i));

    
    Tp = c_pram.period/10.0;
    Hs = c_pram.wave_height/100.0;
    
    Tind  = fn_Tind(1,0.25*Tp,ProbeLocXY(1),WaveType);
    
    t0 = min(Tind(7).time,tm(end));
    
    Tind  = fn_Tind(1,1.75*Tp,ProbeLocXY(1),WaveType);
    t1 = min(Tind(10).time,tm(end));
    
    time = tm(and(tm>t0,tm<t1));
    elev = disp(and(tm>t0,tm<t1));
    figure()
    plot(time,elev)
    [~,~,WTimeD1] = fn_ModifiedMZeroCross([tm,disp]);
    TSD1_IWFC = size(WTimeD1,1);
    
    CalibRef = [CalibRef,TSD1_IWFC];
    tDiff = t1 - t0;
    CalibRefTDiff = [CalibRefTDiff,tDiff];

end
 CalibRefMean = mean(CalibRef);
 CalibRefpersec = CalibRef ./ CalibRefTDiff;
 
MatFile_NM = strcat('Data/Gen/CalibProbe1_',num2str(Tp),'_',num2str(Hs),'.mat');
save(MatFile_NM,'time','elev','Tp','Hs');

%Experiment
for i = 1:size(probesTra,2)
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(39,TNameExp,probesTra(i));

    
    Tp = c_pram.period/10.0;
    Hs = c_pram.wave_height/100.0;
    
    Tind  = fn_Tind(39,0.25*Tp,ProbeLocXY(1),WaveType);
    
    t0 = min(Tind(7).time,tm(end));
    
    Tind  = fn_Tind(39,1.75*Tp,ProbeLocXY(1),WaveType);
    t1 = min(Tind(10).time,tm(end));
    
    time = tm(and(tm>t0,tm<t1));
    elev = disp(and(tm>t0,tm<t1));
    figure()
    plot(time,elev)
    [~,~,WTimeD1] = fn_ModifiedMZeroCross([tm,disp]);
    TSD1_IWFC = size(WTimeD1,1);
    
    ExpTra  = [ExpTra ,TSD1_IWFC];
    tDiff = t1 - t0;
    ExpTraTDiff = [ ExpTraTDiff,tDiff ];

end

ExpTraMean = mean(ExpTra);
ExpTrapersec =  ExpTra ./ ExpTraTDiff;

MatFile_NM = strcat('Data/Gen/ExpProbe11_',num2str(Tp),'_',num2str(Hs),'.mat');
save(MatFile_NM,'time','elev','Tp','Hs');




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

end
