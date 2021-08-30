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

function fn_TransmissionRegAndIrreg(conc,TestName,probes)

if ~exist('PLOT_TYP','var'); PLOT_TYP ='trans_coeff'; end

if ~exist('conc','var'); conc= 39; end ; %79

if ~exist('binnum','var');  binnum=30; end;

if ~exist('ColsReg','var'); ColsReg= '#058000'; end ;
if ~exist('Symbols','var'); Symbols= {'x','d','^','s'}; end ;
if ~exist('Cols','var');Cols= {'#ff0000','#0bff01 ','#0487f9','#9701ff'}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=2e2; end
if ~exist('model_pers','var'); model_pers=0.5:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=1; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;


%Load wikiwaves
WikiWavesPRSATransmission = load('./Data/Gen/WikiWavesElasticTransmission');


Description = ['Concentration ' num2str(conc) '% '];



 th_res=100;
SURGE=0;
terms_grn=100;
extra_pts=[];     
rigid = 4; 
N=800; %determines the number of points along beam (i.e. floating plate)
n=400; %number of eigenvalues
nDS = n;
Param = ParamDef_Oceanide(rigid); 
Param = ModParam_def(Param,1,n,extra_pts,terms_grn,th_res); 

Param.rho = WikiWavesPRSATransmission.PhysVars.rho;
Param.E = WikiWavesPRSATransmission.PhysVars.E;
Param.thickness = WikiWavesPRSATransmission.PhysVars.thickness;
Param.bed = WikiWavesPRSATransmission.PhysVars.h;
Param.floe_diam = 2*WikiWavesPRSATransmission.PhysVars.L;

Param.draft = Param.thickness.*Param.rho./Param.rho_0; 

% Flexural rigidities (in MPa\,m^3)
Param.D = Param.E.*Param.thickness.^3./(12*(1-Param.nu.^2));

% Scaled rigidity 
Param.beta = Param.D./(Param.g.*Param.rho_0);

TpS = 0.5:0.02:2;

for i = 1:length(TpS)
    out = fn_ElasticRaft2d('freq',1/TpS(i),Param,'disk',0,0,1,1,0,'r',[0]);
    T(i) = abs(out(1).value);
    R(i) = abs(out(2).value);
end


figure('DefaultAxesFontSize',18);
hold on;
clear PlotNames;
PlotNames = {};

%models
plot(WikiWavesPRSATransmission.Tps,WikiWavesPRSATransmission.T_GFA,'-b','LineWidth',2);
pname = 'Greens Function Wiki Waves';
PlotNames{end +1} = pname;

plot(WikiWavesPRSATransmission.Tps,WikiWavesPRSATransmission.T_EMA,'--k','LineWidth',2);
pname = 'Eigenmode Match Wiki Waves';
PlotNames{end +1} = pname;

plot(TpS,T,'-r','LineWidth',2);
pname = 'Luke Variational';
PlotNames{end +1} = pname;

xlabel('Period (s)');
ylabel('|T|) ');
title(['|T| ']);
legend(PlotNames);
xlim([0 , 2 ]);
ylim([0 , 1 ]);



 return

