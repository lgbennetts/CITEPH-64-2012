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

function fn_TransmissionNewAvg(conc,TestName,probes)

% %Compare calibration to disk


if ~exist('conc','var') conc={0}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {'14'}; end ; 
if ~exist('wbin','var'); wbin =3; end ;  

% if ~exist('Cols','var'); Cols= {'xr', '^b','sk','dg'}; end ;
if ~exist('Cols','var'); Cols= {{'#ff0000','#680303'},{'#0bff01 ','#058000'},{'#0487f9','#00427c'},{'#9701ff','#470179'}}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

if ~exist('PerNum','var');  PerNum=7; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end
if ~exist('OutputStr','var');  OutputStr = 'TargetSpectra'; end

if ~exist('t0description','var'); t0description ='waves reach x'; end ; 
if ~exist('t1description','var'); t1description ='final waves reach x'; end ; 


%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;

for j = 1: length(TestName)
    
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}, probes{1});
    Tp = c_pram.period/10.0;
    Hs = c_pram.wave_height/100.0;
    
    dt = tm(2) - tm(1); 
     
    Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),'Irregular');
    t0index = find(strcmp({Tind.description},t0description));
    t0 = Tind(t0index).time; %+ T_windowSec/2;
    t1index = find(strcmp({Tind.description},t1description));
    t1 = Tind(t1index).time;% - T_windowSec/2 ;
    
    TI = 50;
    %
     figure('DefaultAxesFontSize',18);
     plot(tm(1:TI:end),disp(1:TI:end),'-k', 'DisplayName','Time Series Calibration Probe 11');
     hold on;
      plot([t0,t0],[-0.02 0.02],'--b', 'DisplayName','First Waves With \tau=0.8s Reach Probe');
      plot([t1,t1],[-0.02 0.02],'--b', 'DisplayName','Final Waves With \tau=0.8s Reach Probe');
     axis([0 550 -0.02 0.02 ])
     title(['Elevation Time Series for Irregular Experiment with \tau = ', num2str(Tp), ' Hs = ', num2str(Hs)]);
     xlabel('time (s)')
     ylabel('elevation (m)')
     legend();
%      matlab2tikz('ElevationTimeSeriesEx.tex'); 
     
     
     %Energy in 
     
    eval(TpersStr);
    T_windowSec = round(Tpers*Tp);
    SamplingRate = floor(1./dt);
    T_window = SamplingRate*T_windowSec;
     
    [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
    period = 1./ff;
    FourierWindow_Sec = FourierWindow/SamplingRate;
    
    
    wJS = 2*pi/Tp;
    JSspec  = jonswap(2*pi*ff,'wp',wJS,'Hs',Hs);
    
    jj0=find(tm_vec>t0,1);
    jj1=find(tm_vec<t1,1,'last');
    [~,Tpj]=min(abs(period-Tp));
    SpecFilt =(1./(period(Tpj + wbin) - period(Tpj - wbin))).* trapz(period(Tpj + [-wbin:wbin]),Sn_mat(Tpj + [-wbin:wbin],:));

    figure('DefaultAxesFontSize',18);
    plot(tm_vec,SpecFilt,'-k', 'DisplayName','Energy At \tau = 0.8');
    hold on;
    plot([0,500],[mean(SpecFilt(jj0:jj1)),mean(SpecFilt(jj0:jj1))],'--g', 'DisplayName','Mean Over Window');
    plot([0,500],[JSspec(Tpj),JSspec(Tpj)],'-r', 'DisplayName','Jonswap Spectrum Energy at \tau=0.8');

    
    plot([t0,t0],[0 2e-4],'--b', 'DisplayName','First Waves With \tau=0.8s Reach Probe');
    plot([t1,t1],[0 2e-4],'--b', 'DisplayName','Final Waves With \tau=0.8s Reach Probe');
    axis([0 500 0 2e-4 ])
    title(['Energy Over Time \tau = ', num2str(Tp), ' Hs = ', num2str(Hs)]);
    xlabel('Mean Time (s)')
    ylabel('Energy (m^2 s^-2)')
    legend();
    
%     matlab2tikz('EnergyAtPeriod.tex'); 
    
    
    TpS = 0.25*Tp:0.05:2*Tp;
    
    PerWholeCalib = [];
    SpectraWholeCalib = [];
    PerWholeConc = [];
    SpectraWholeConc = [];
    for ji = 1: length(TpS)
        OutputStr = 'TargetSpectra';
%         outCalib = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{1},TestName{j}{1},probes{1},OutputStr);
        outCalib = fn_SpectraAverageTgT(TpS(ji),conc{1},TestName{j},WaveType,11,OutputStr);
        PerFilt = outCalib(1).value;
        SpectraFilt = outCalib(2).value;
        
        PerWholeCalib = [PerFilt;PerWholeCalib];
        SpectraWholeCalib = [SpectraFilt;SpectraWholeCalib];      
    end

    ffWholeConc = 1./PerWholeConc;
    
    wJS = 2*pi/Tp;
    ffJS = 0:0.001:10*Tp;
    JSspec  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
    
    %Figure
     figure();
     plot(PerWholeCalib,SpectraWholeCalib,'-b', 'DisplayName','Calibration');
     hold on;
     plot(PerWholeCalib,SpectraWholeCalib,'-b', 'DisplayName','Calibration');
     plot(1./ffJS,JSspec,'--k' , 'DisplayName','JonSwap');
     axis([0.1*Tp 2*Tp 0 7e-5 ])
     xlabel('T(s)')
     ylabel('S(T)')
%      matlab2tikz('SpectraSingleProbe.tex');
     
     
     %Multiple Probes
         
    PerWholeCalib = [];
    SpectraWholeCalib = [];
    PerWholeConc = [];
    SpectraWholeConc = [];
    for ji = 1: length(TpS)
        OutputStr = 'TargetSpectra';
%         outCalib = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{1},TestName{j}{1},probes{1},OutputStr);
        outCalib = fn_SpectraAverageTgT(TpS(ji),conc{1},TestName{j},WaveType,11:20,OutputStr);
        PerFilt = outCalib(1).value;
        SpectraFilt = outCalib(2).value;
        
        PerWholeCalib = [PerFilt;PerWholeCalib];
        SpectraWholeCalib = [SpectraFilt;SpectraWholeCalib];      
    end

    ffWholeConc = 1./PerWholeConc;
    
%      figure();
%      plot(PerWholeCalib,SpectraWholeCalib,'-b', 'DisplayName','Calibration');
%      hold on;
%      plot(PerWholeCalib,SpectraWholeCalib,'-b', 'DisplayName','Calibration');
%      plot(1./ffJS,JSspec,'--k' , 'DisplayName','JonSwap');
%      axis([0.1*Tp 2*Tp 0 7e-5 ])
%      xlabel('T(s)')
%      ylabel('S(T)')
     matlab2tikz('SpectraMultipleProbe.tex'); 

   
    
    
    j = j +1;
end

%Save Matrices
% TpTarg = TpA;
% HsTarg = HsA;
% F = ffFilt;
% WL =  wlExp;
% Trans = Tran;
% Ain = ATWCalibFilt;
% 
% MatFile_NM = strcat('Data/Gen/TransEnergy',int2str(conc{2}),WaveType(1:3));
% save(MatFile_NM,'TpTarg','HsTarg','F','WL','Trans','Ain');


function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

