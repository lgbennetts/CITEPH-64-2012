%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

clear all
close all
% clc

%currently configured to run multiple tests via these opening for loops.

%select wave period (s)

waveperiod = 0.95


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

Es = 0.5:0.5:12;

for i = 1:length(Es)
 Es(i)
Param.E = 10^Es(i);
% Flexural rigidities (in MPa\,m^3)
Param.D = Param.E.*Param.thickness.^3./(12*(1-Param.nu.^2));

% Scaled rigidity 
Param.beta = Param.D./(Param.g.*Param.rho_0);
xbeam = [0]
out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',SURGE,0,1,1,0,'r',xbeam);
T_EM_abs(i)  = abs(out(1).value);
R_EM_abs(i)  = abs(out(2).value);
EtaOut = out(3).value;
uSG = out(4).value;



 
end


figure('DefaultAxesFontSize',18);
hold on;
clear PlotNames;
PlotNames = {};

%models
plot(Es,R_EM_abs,'.b','MarkerSize',12);
pname = '|R|  Luke';
PlotNames{end +1} = pname;

plot(Es,T_EM_abs,'.r','MarkerSize',12);
pname = '|T|  Luke';
PlotNames{end +1} = pname;


plot([1,10],[abs(0.3570),abs(0.3570)],'--b','LineWidth',3);
pname = '|R|  Lucas';
PlotNames{end +1} = pname;

plot([1,10],[abs( 0.9341),abs( 0.9341)],'--r','LineWidth',3);
pname = '|T|  Lucas';
PlotNames{end +1} = pname;

xlabel('Power of E');
ylabel('|Value|');
title(['Rigid Convergence PRSA Properties' ]);
legend(PlotNames);








