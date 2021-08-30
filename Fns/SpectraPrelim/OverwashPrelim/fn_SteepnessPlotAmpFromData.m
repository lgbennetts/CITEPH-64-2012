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

function fn_SteepnessPlot()

if ~exist('conc','var') conc =39; end %39,79 or empty (1)
if ~exist('TestNames','var') TestNames ={'2','2a','2b','7'}; end 
if ~exist('TestNames','var') probes =11:20; end 

if ~exist('TgTp','var') TgTp =0.95; end 
if ~exist('WaveType','var') WaveType = 'Regular'; end 


close all;
%Calibration
outRef = fn_SpectraAverageTgT(TgTp,0,'4',1:10,'TargetSpectra');
CalibRef = sqrt(2*outRef(2).value);
outTrans = fn_SpectraAverageTgT(TgTp,0,'4',11:20,'TargetSpectra');
CalibTran = sqrt(2*outTrans(2).value);

CalibT = [1,1,1,2]*CalibTran;

%Ok Lets Loop through TestNames
for j = 1 : length(TestNames)
    outRef = fn_SpectraAverageTgT(TgTp,conc,TestNames{j},1:10,'TargetSpectra');
    ExpRef(j) = sqrt(2*outRef(2).value);   
    outTrans = fn_SpectraAverageTgT(TgTp,conc,TestNames{j},11:20,'TargetSpectra');
    ExpTrans(j) = sqrt(2*outTrans(2).value);
    
    if strcmp(TestNames{j},'7')
       TCE(j) = ExpTrans(j) / (2*CalibTran); 
    else
       TCE(j) = ExpTrans(j) / (CalibTran);
    end
    
    ToRp1(j) = ExpTrans(j) / ExpRef(j);
    
    TS(j) = 2*ToRp1(j) / (1 +ToRp1(j));
    
    
end



% %Find incoming and outgoing wave steepness
% 
if conc == 39
   Mat = load('./Data/Trans38.mat');
else
   Mat = load('./Data/Trans77.mat'); 
end


IndexMatchTp = find([Mat.T] == TgTp) ; 

HsList = [Mat.H];
HsList = HsList(IndexMatchTp);
AsList = HsList/1000.0;

ETransList = [Mat.trans];
ETransList = ETransList(IndexMatchTp);

ATrans = sqrt(ETransList);


%Overwash - Skene output

OWData = load('Tp0.95SingleFloeOW.mat');

%k
k = 4.45906;

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=[0.95]; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

%Models - Transmission Single Floe
th_res=100;
SURGE=1;
terms_grn=100;
extra_pts=[];

if ~exist('Param','var'); Param = ParamDef_Oceanide(4); 
Param = ModParam_def(Param,1,1e2,extra_pts,terms_grn,th_res); end

%Linear Theories
if conc==39
dum_c = pi*((Param.floe_diam/2)^2)/2;
elseif conc==79
dum_c = pi*((Param.floe_diam/2)^2);
end
EDN_Real = (dum_c*Param.MIZ_length)/Param.floe_diam;
EDN_2D = (dum_c*Param.MIZ_length)/Param.floe_diam/2;

out = fn_ElasticRaft2d('freq',1./TgTp,Param,'transmitted energy',SURGE,0,1,1,0,'r');
Trans = out.value(1);
% TransMultiple = (exp(dum_c/100.0*log(Trans)*Param.MIZ_length/Param.floe_diam /2));   

TwoDEMM = Main_AttnModels(model_pers,dum_c*100,'2d EMM',0,Vert_Modes,DO_FDSP,0);
TwoDEMM = (TwoDEMM.value);

x = 0: 0.1:1;

%OW
OWtransT = OWData.kT ./OWData.kI;
OWtransT(isnan(OWtransT)) = OWtransT(end);
% 
OWtransTManyDisks = OWtransT.^EDN_Real.*OWData.kI;

% OW
OWtransTN = OWData.kTnew ./ OWData.kI;
OWtransTN(isnan(OWtransTN)) = OWtransT(end);
% 
InterExtrap = [0:0.01:2,OWData.kI(2:end)];
InterExtrap = sort(InterExtrap);
TOW1 = interp1(OWData.kI(2:end),OWData.kTnew(2:end),InterExtrap,'spline','extrap');
TOW2 = interp1(OWData.kI(2:end),OWData.kTnew(2:end),TOW1,'spline','extrap');

% one overwash then, regular transmission
MultipleDiskTrans =OWtransTN.^EDN_Real;
OWtransTNManyDisks = MultipleDiskTrans.*OWData.kI;

figure();
set(gca,'FontSize',18) 
% plot(k*AsList,k*(ATrans).*AsList,'xb','MarkerSize',16, 'DisplayName', 'Data');
plot(k*CalibT,k*ExpTrans,'xb','MarkerSize',16, 'DisplayName', 'Data');
hold on;
plot(k*(ExpRef)./(2 - TS),k.*ExpTrans,'xr','MarkerSize',16, 'DisplayName', 'Data 1/(r+1)');
plot(k*ExpRef,k*ExpTrans,'xk','MarkerSize',16, 'DisplayName', 'Data Ref/Trans');

plot(k*AsList,k*(ATrans).*AsList,'ob','MarkerSize',16, 'DisplayName', 'Paper Analysis');
% plot(OWData.kI,OWData.kT,'-r', 'DisplayName', 'Single Floe Model','LineWidth',2);
plot(OWData.kI,OWtransTManyDisks,'-r', 'DisplayName', 'Linear Many Floe Model','LineWidth',2);
plot(InterExtrap,TOW1,'-k', 'DisplayName', 'Extrap Many Floe Model Overwash - First Iteration (One Floe)','LineWidth',2);
plot(OWData.kI,OWData.kTnew,'.k', 'DisplayName', 'OW Floe Data','MarkerSize',12);

plot(InterExtrap,TOW2,'-g', 'DisplayName', 'Extrap Many Floe Model Overwash - Second Iteration (Second Floe)','LineWidth',2);% plot(OWData.kI,OWtransTN.^effectivedisknumber.*OWData.kI,'--k', 'DisplayName', 'Extrap Many Floe Model Overwash -  Power');
plot(OWData.kI,OWtransTNManyDisks,'--g', 'DisplayName', 'Many Floe Model Overwash (simple Power)','LineWidth',2);

title(['Tp = ', num2str(TgTp), ' cL/a: ', num2str(EDN_Real)]);
xlabel('Incoming Wave Steepeness');
ylabel('Transmitted Wave Steepeness');
xlim([0,0.5])
legend();


%normalise to do a better iteration job
return



