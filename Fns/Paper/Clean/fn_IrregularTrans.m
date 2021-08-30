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

function fn_IrregularTrans(conc,TestName,probes)

% %Compare calibration to disk
if ~exist('conc','var') conc={0,39}; end %39,79 or empty (1)
if ~exist('probes_tra','var'); probes_tra= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('probes_ref','var'); probes_ref= {1:10,1:10} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  
%if ~exist('TestName','var'); TestName= {{'14','9'},{'15','10'},{'16','8'}}; end ;  

% if ~exist('conc','var') conc={0,79}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% % if ~exist('TestName','var'); TestName= {{'14','19'}}; end ; 
% if ~exist('TestName','var'); TestName= {{'14','19'}, {'15','20'},{'16','21'},{'17','22'}}; end ; 


% if ~exist('Cols','var'); Cols= {'xr', '^b','sk','dg'}; end ;
if ~exist('Cols','var'); Cols= {{'#ff0000','#680303'},{'#0bff01 ','#058000'},{'#0487f9','#00427c'},{'#9701ff','#470179'}}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

if ~exist('PerNum','var');  PerNum=30; end %28; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end
if ~exist('OutputStr','var');  OutputStr = 'TargetSpectra'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;

for j = 1: length(TestName)
    
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, 1);
    Tp = c_pram.period/10.0;
    Hs = c_pram.wave_height/100.0;
    
%     TpS = linspace(0.75*Tp,Tp,PerNum);
%     TpS = [TpS(1:end-1),linspace(Tp,1.5*Tp,PerNum)];
    
    TpS = 0.25*Tp :Tp /PerNum: 1.75*Tp;
    
    PerWholeCalib = [];
    SpectraWholeCalib = [];
    PerWholeConc = [];
    SpectraWholeConc = [];
    
    PerWholeCalib_ref  = [];
    SpectraWholeCalib_ref  = [];
    PerWholeConc_ref  = [];
    SpectraWholeConc_ref  = [];
    for ji = 1: length(TpS)
        OutputStr = 'TargetSpectra';
%         outCalib = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{1},TestName{j}{1},probes{1},OutputStr);
        outCalib = fn_SpectraAverageTgT(TpS(ji),conc{1},TestName{j}{1},WaveType,probes_tra{1},OutputStr);
        PerFilt = outCalib(1).value;
        SpectraFilt = outCalib(2).value;
        
        PerWholeCalib = [PerFilt;PerWholeCalib];
        SpectraWholeCalib = [SpectraFilt;SpectraWholeCalib];
        
%         outConc = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{2},TestName{j}{2},probes{2},OutputStr);
        outConc = fn_SpectraAverageTgT(TpS(ji),conc{2},TestName{j}{2},WaveType,probes_tra{2},OutputStr);
        PerFilt = outConc(1).value;
        SpectraFilt = outConc(2).value;
        
        PerWholeConc = [PerFilt;PerWholeConc];
        SpectraWholeConc = [SpectraFilt;SpectraWholeConc];  
        
        
        outCalib_ref = fn_SpectraAverageTgT(TpS(ji),conc{1},TestName{j}{1},WaveType,probes_ref{1},OutputStr);
        PerFilt_ref  = outCalib_ref(1).value;
        SpectraFilt_ref  = outCalib_ref(2).value;
        
        PerWholeCalib_ref  = [PerFilt_ref ;PerWholeCalib_ref ];
        SpectraWholeCalib_ref  = [SpectraFilt_ref ;SpectraWholeCalib_ref ];
        
%         outConc = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{2},TestName{j}{2},probes{2},OutputStr);
        outConc_ref = fn_SpectraAverageTgT(TpS(ji),conc{2},TestName{j}{2},WaveType,probes_ref{2},OutputStr);
        PerFilt_ref = outConc_ref(1).value;
        SpectraFilt_ref = outConc_ref(2).value;
        
        PerWholeConc_ref = [PerFilt_ref;PerWholeConc_ref];
        SpectraWholeConc_ref = [SpectraFilt_ref;SpectraWholeConc_ref];    
    end

    ffWholeConc = 1./PerWholeConc;
    
    wJS = 2*pi/Tp;
    ffJS = 0:0.01:10;
    JSspec{j}  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
    
    %Figure - Spectra
%      figure();
%      plot(PerWholeCalib,SpectraWholeCalib,'-b', 'DisplayName','Calibration');
%      hold on;
%      plot(PerWholeConc,SpectraWholeConc,'-r', 'DisplayName','Experiment');
%      plot(1./ffJS,JSspec{j},'--k' , 'DisplayName','JonSwap');
%      title(['Tp = ', num2str(Tp), ' Hs = ', num2str(Hs)]);
%      xlim([0 2*Tp])
%      xlabel('T(s)')
%      ylabel('Spectra')
    
   
    %experimental
    TpA{j} = Tp;
    HsA{j} = Hs;
    ffA{j} = ffWholeConc;
    ffpA{j} =ffWholeConc*Tp;
    
   
    ATWCalib_Tra{j} = SpectraWholeCalib;
    ATWConc_Tra{j} = SpectraWholeConc;
    
    ATWCalib_Ref{j} = SpectraWholeCalib_ref;
    ATWConc_Ref{j} = SpectraWholeConc_ref;

    Tran{j} = (ATWConc_Tra{j}./ATWCalib_Tra{j});

    clear wlc;
    for i = 1: length(ffA{j})
        ki = dispersion_free_surface((2*pi*ffA{j}(i))^2/9.81,0,3.1);
        wlc(i,1) = 2*pi/imag(ki);
    end
    
    wlExp{j} = wlc;
    
    j = j +1;
end

%Save Matrices
TpTarg = TpA;
HsTarg = HsA;
F = ffA;
WL =  wlExp;
Trans = Tran;
% Ein = ATWCalibFilt;
CalibTra = ATWCalib_Tra;
CalibRef = ATWCalib_Ref;
ExpTra = ATWConc_Tra;
ExpRef = ATWConc_Ref;


MatFile_NM = strcat('Data/Gen/NewTransEnergy',int2str(conc{2}),WaveType(1:3));
save(MatFile_NM,'TpTarg','HsTarg','F','WL','Trans','CalibTra','CalibRef','ExpTra','ExpRef');



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

