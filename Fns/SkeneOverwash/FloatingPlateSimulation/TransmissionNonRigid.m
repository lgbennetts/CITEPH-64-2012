%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

clear all
close all
% clc

%currently configured to run multiple tests via these opening for loops.

%select wave period (s)




th_res=100;
SURGE=0;
terms_grn=100;
extra_pts=[];     
rigid = 4; 
N=100; %determines the number of points along beam (i.e. floating plate)
n=100; %number of eigenvalues
nDS = n;
if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid); 
    Param = ModParam_def(Param,1,n,extra_pts,terms_grn,th_res); end


Tps = 0.5:0.02:2;

%Load wikiwaves
WikiWavesPRSATransmission = load('./Data/Gen/WikiWavesElasticTransmission');

for i = 1:length(Tps)
% Scaled rigidity 
Param.beta = Param.D./(Param.g.*Param.rho_0);
xbeam = [0];
out = fn_ElasticRaft2d('freq',1/Tps(i),Param,'disk',SURGE,0,1,1,0,'r',xbeam);
T_EM_abs(i)  = abs(out(1).value)^2;
R_EM_abs(i)  = abs(out(2).value);
EtaOut = out(3).value;
uSG = out(4).value;



 
end



figure('DefaultAxesFontSize',18);
hold on;
clear PlotNames;
PlotNames = {};

%models
plot(Tps,T_EM_abs,'.k','MarkerSize',12);
pname = '|T|  Luke';
PlotNames{end +1} = pname;


plot(WikiWavesPRSATransmission.Tps,WikiWavesPRSATransmission.T_GF_Energy,'--b','LineWidth',2);
pname = 'Greens Function Wiki Waves';
PlotNames{end +1} = pname;

plot(WikiWavesPRSATransmission.Tps,WikiWavesPRSATransmission.T_EM_Energy,'--r','LineWidth',2);
pname = 'Eigen Wiki Waves';
PlotNames{end +1} = pname;

xlabel('period');
ylabel('transmission');
title(['Elastic Convergence PRSA Properties' ]);
legend(PlotNames);








