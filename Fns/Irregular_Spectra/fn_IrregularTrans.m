% function fn_IrregularTrans
%
% DESCRIPTION: Generate the transmissions matrices
%
% INPUTS:
%
% OUTPUTS:
%     MatFile_NM - saved as a matrix file
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by

function fn_IrregularTrans(conc)

% %Compare calibration to disk
if ~exist('probes_tra','var'); probes_tra= {11:20,11:20} ;end; 
if ~exist('probes_ref','var'); probes_ref= {1:10,1:10} ;end;  
if ~exist('conc','var'); conc= 39 ;end; 

if conc == 39
    conc_tests={0,39}; 
    TestName={{'14','9'}};
%     TestName= {{'14','9'},{'15','10'},{'16','8'}};
elseif conc == 79
    conc_tests={0,79}; 
    TestName= {{'14','19'}, {'15','20'},{'16','21'},{'17','22'}};
end

if ~exist('Cols','var'); Cols= {{'#ff0000','#680303'},{'#0bff01 ','#058000'},{'#0487f9','#00427c'},{'#9701ff','#470179'}}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.01:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

if ~exist('PerNum','var');  PerNum=100; end %28; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end
if ~exist('OutputStr','var');  OutputStr = 'TargetSpectra'; end



for j = 1: length(TestName)
    
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc_tests{1},TestName{j}{1}, 1);
    Tp = c_pram.period/10.0;
    Hs = c_pram.wave_height/100.0;
    
    TpS = linspace(0.1,4*Tp,PerNum);
%     TpS = [TpS(1:end-1),linspace(Tp,1.5*Tp,PerNum)];
    
%     TpS = 0.25*Tp :Tp /PerNum: 1.75*Tp;

%     TpS = linspace(0.1,Tp,PerNum);
    
    PerWholeCalib = [];
    SpectraWholeCalib_A = [];
    SpectraWholeCalib_S = [];
    PerWholeConc = [];
    SpectraWholeConc_A = [];
    SpectraWholeConc_S = [];
    
    PerWholeCalib_ref  = [];
    SpectraWholeCalib_A_ref  = [];
    SpectraWholeCalib_S_ref  = [];
    PerWholeConc_ref  = [];
    SpectraWholeConc_A_ref  = [];
     SpectraWholeConc_S_ref  = [];
    for ji = 1: length(TpS)
        OutputStr = 'TargetSpectra';
%         outCalib = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{1},TestName{j}{1},probes{1},OutputStr);
        outCalib = fn_SpectraAverageAtTgPer(TpS(ji),conc_tests{1},TestName{j}{1},WaveType,probes_tra{1},OutputStr);
        PerFilt = outCalib(1).value;
        SpectraFilt_A = outCalib(2).value;
        SpectraFilt_S = outCalib(3).value;
        
        PerWholeCalib = [PerFilt;PerWholeCalib];
        SpectraWholeCalib_A = [SpectraFilt_A;SpectraWholeCalib_A];
        SpectraWholeCalib_S = [SpectraFilt_S;SpectraWholeCalib_S];
        
%         outConc = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{2},TestName{j}{2},probes{2},OutputStr);
        outConc = fn_SpectraAverageAtTgPer(TpS(ji),conc_tests{2},TestName{j}{2},WaveType,probes_tra{2},OutputStr);
        PerFilt = outConc(1).value;
        SpectraFilt_A = outConc(2).value;
        SpectraFilt_S = outConc(3).value;
        
        PerWholeConc = [PerFilt;PerWholeConc];
        SpectraWholeConc_A = [SpectraFilt_A;SpectraWholeConc_A];  
        SpectraWholeConc_S = [SpectraFilt_S;SpectraWholeConc_S]; 
        
        
        outCalib_ref = fn_SpectraAverageAtTgPer(TpS(ji),conc_tests{1},TestName{j}{1},WaveType,probes_ref{1},OutputStr);
        PerFilt_ref  = outCalib_ref(1).value;
        SpectraFilt_A_ref  = outCalib_ref(2).value;
        SpectraFilt_S_ref  = outCalib_ref(3).value;
        
        PerWholeCalib_ref  = [PerFilt_ref ;PerWholeCalib_ref ];
        SpectraWholeCalib_A_ref  = [SpectraFilt_A_ref ;SpectraWholeCalib_A_ref ];
        SpectraWholeCalib_S_ref  = [SpectraFilt_S_ref ;SpectraWholeCalib_S_ref ];
        
%         outConc = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{2},TestName{j}{2},probes{2},OutputStr);
        outConc_ref = fn_SpectraAverageAtTgPer(TpS(ji),conc_tests{2},TestName{j}{2},WaveType,probes_ref{2},OutputStr);
        PerFilt_ref = outConc_ref(1).value;
        SpectraFilt_A_ref = outConc_ref(2).value;
        SpectraFilt_S_ref = outConc_ref(3).value;
        
        PerWholeConc_ref = [PerFilt_ref;PerWholeConc_ref];
        SpectraWholeConc_A_ref = [SpectraFilt_A_ref;SpectraWholeConc_A_ref]; 
        SpectraWholeConc_S_ref = [SpectraFilt_S_ref;SpectraWholeConc_S_ref];
    end

    ffWholeConc = 1./PerWholeConc;
    
    wJS = 2*pi/Tp;
    ffJS = 0:0.01:10;
    JSspec{j}  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
    
    %Figure - Spectra
     figure();
     errorbar(PerWholeCalib,SpectraWholeCalib_A,SpectraWholeCalib_S,'-b', 'DisplayName','Calibration');
     hold on;
     errorbar(PerWholeConc,SpectraWholeConc_A,2*SpectraWholeConc_S,'-r', 'DisplayName','Experiment');
     plot(1./ffJS,JSspec{j},'--k' , 'DisplayName','JonSwap');
     title(['Tp = ', num2str(Tp), ' Hs = ', num2str(Hs)]);
     xlim([0 2*Tp])
     xlabel('T(s)')
     ylabel('Spectra')
    
   
    %experimental
    TpA{j} = Tp;
    HsA{j} = Hs;
    ffA{j} = ffWholeConc;
    ffpA{j} =ffWholeConc*Tp;
    
   
    ATWCalib_A_Tra{j} = SpectraWholeCalib_A;
    ATWCalib_S_Tra{j} = SpectraWholeCalib_S;
    ATWConc_A_Tra{j} = SpectraWholeConc_A;
    ATWConc_S_Tra{j} = SpectraWholeConc_S;
    
    ATWCalib_A_Ref{j} = SpectraWholeCalib_A_ref;
    ATWCalib_S_Ref{j} = SpectraWholeCalib_S_ref;
    ATWConc_A_Ref{j} = SpectraWholeConc_A_ref;
    ATWConc_S_Ref{j} = SpectraWholeConc_S_ref;

    Tran{j} = (ATWConc_A_Tra{j}./ATWCalib_A_Tra{j});
    TranSTD{j} = (ATWConc_S_Tra{j}./ATWCalib_A_Tra{j});

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
TransSTD = TranSTD;
% Ein = ATWCalibFilt;
CalibTra_A = ATWCalib_A_Tra;
CalibRef_A = ATWCalib_A_Ref;
ExpTra_A = ATWConc_A_Tra;
ExpRef_A = ATWConc_A_Ref;
CalibTra_S = ATWCalib_S_Tra;
CalibRef_S = ATWCalib_S_Ref;
ExpTra_S = ATWConc_S_Tra;
ExpRef_S = ATWConc_S_Ref;


MatFile_NM = strcat('Data/Gen/A_NewTransAmp',int2str(conc_tests{2}),WaveType(1:3));
save(MatFile_NM,'TpTarg','HsTarg','F','WL','Trans','CalibTra_A','CalibRef_A','ExpTra_A','ExpRef_A','TransSTD','CalibTra_S','CalibRef_S','ExpTra_S','ExpRef_S');



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

