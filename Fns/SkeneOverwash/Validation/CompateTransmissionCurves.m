%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

clear all
close all
% clc

if ~exist('ColsReg','var'); ColsReg= '#058000'; end ;
if ~exist('Symbols','var'); Symbols= {'x','d','^','s'}; end ;
if ~exist('Cols','var');Cols= {'#ff0000','#0bff01 ','#0487f9','#9701ff'}; end ;

%currently configured to run multiple tests via these opening for loops.
WavePer = 0.5:0.02:2;
TDS = zeros(size(WavePer));
TLB = zeros(size(WavePer));
TLBNS = zeros(size(WavePer));
WL = zeros(size(WavePer));
TLBLFNS = zeros(size(WavePer));

th_res=100;
SURGE=1;
terms_grn=100;
extra_pts=[];     
rigid = 4; 
N=100; %determines the number of points along beam (i.e. floating plate)
Nd = N;
n=100; %number of eigenvalues
nd = n;
if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid); 
    Param = ModParam_def(Param,1,n,extra_pts,terms_grn,th_res); end
%  
    
E=Param.E; %plate elasticity Pa
rho_i=Param.rho; %plate density kg/m^3
nu=Param.nu; %poisson's ratio

thickness=Param.thickness; %plate thinkcness m
rho_w=Param.rho_0; %water density kg/m^3
g=9.81; %gravitational accel m/s^2
H=Param.bed; %water depth m
L=Param.floe_diam; %beam length m

draughtratio=rho_i/rho_w ;%archemedies submergence
plateabove=thickness*(1-draughtratio); %portion of plate above water
platebelow=thickness*draughtratio; %portion of plate below water

%select wave period (s)
for i = 1:length(WavePer)
    %select incident wave steepness
    waveperiod=WavePer(i);

    freq=1/waveperiod; %wave frequency Hz
    omega=2*pi*freq; %angular frequency rad/s

    alpha=omega^2/g; %scaled frequency
    k=dispersion_free_surface(alpha,0,H); %wavenumber
    kr = imag(k);

    %configuring some nondimensional variables relevant to LPT solution
%      gamma=rho_i*thickness/rho_w;
%      beta=E*thickness^3/(12*(1-nu^2)*rho_w/g);

%     D = E*thickness^3/(12*(1-nu^2));
%     CL = (D / (rho_w*omega^2))^(1/5);
%     alpha = CL*omega^2/g;
%     beta = D / (rho_w*g*CL^4);
%     gamma = rho_i*thickness/(rho_w*CL);

    D = E*thickness^2/(12*rho_i*(1-nu^2));
    CL = L;
    alpha = CL*omega^2/g;
    beta = (rho_i*thickness*D) / (rho_w*g*CL^4);
    gamma = rho_i*thickness/(rho_w*CL);
    
    %Get linear potential theory solution
    %displacement is plate's frequency domain displacement
    %R and T are the reflected and transmitted complex amplitudes that are
    %scaled
    [xi,displacement,potential,R,T,f] = elastic_plate_modes(alpha,beta,gamma,H,L/2,nd,N);
    xbeam=-L/2:L/(length(displacement)-1):L/2;

    out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',1,0,0,1,0,'r',xbeam);
    TN = (out(1).value);
    RN = (out(2).value);
    
%     out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',0,0,0,1,0,'r',xbeam);
%     TNNS = (out(1).value);
%     RNNS = (out(2).value);
%     
%     out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'transmitted energy',0,1,0,1,0,'r',xbeam);
%     TNLF = (out(1).value);
%     RNLF = 1- TNLF;

    TDS(i) = abs(T)^2;
    TLB(i) = abs(TN)^2;
%     TLBNS(i) = abs(TNNS)^2;
%     TLBLFNS(i) = TNLF;
%     RDS(i) = abs(R)^2;
%     RLB(i) = abs(RN)^2;
%     RLBNS(i) = abs(RNNS)^2;
    WL(i) = 2*pi/kr;
    


end

%basic
figure()
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};
plot(WL,TDS,'-r')
PlotNames{end +1} = 'Wiki Waves';
hold on;
plot(WL,TLB,'-k')
PlotNames{end +1} = 'Lukes';
% plot(WL,TLBNS,'-b')
% PlotNames{end +1} = 'Lukes No Surge';
xlabel('Wavelength (m)')
ylabel('Transmission Coefficient (Energy)')
title('Single Floe Transmission')
legend(PlotNames);


c = pi*(0.495^2)/2;
MultTDS =  TDS.^(c*Param.MIZ_length/(Param.floe_diam));
MultTLB = TLB.^(c*Param.MIZ_length/(Param.floe_diam));
MultTLBNS = TLBNS.^(c*Param.MIZ_length/(Param.floe_diam));
MultTLBLFNS = TLBLFNS.^(c*Param.MIZ_length/(Param.floe_diam));

Boltz = Main_AttnModels(WavePer,100*c,'Boltzmann steady',0,n,0,0);
Boltz = (Boltz.value).^2;

% TransDataIrreg = load('./Data/Gen/NewTrans39Irr');
TransDataReg = load('./Data/Trans38');
for i = 1: length(TransDataReg.T)
    ki = dispersion_free_surface((2*pi/TransDataReg.T(i))^2/Param.g,0,H);
    WLReg(i) = 2*pi/imag(ki);
end

TransDataRegNew = load('./Data/Gen/TransEnergy39Reg');
for i = 1: length(TransDataRegNew.TpA)
    ki = dispersion_free_surface((2*pi/TransDataRegNew.TpA(i))^2/Param.g,0,H);
    WLRegNew(i) = 2*pi/imag(ki);
end

%irregular
TransDataIrreg = load('./Data/Gen/TransEnergy39Irr');
for i = 1:length(TransDataIrreg.Trans)
    
    for j = 1:length(TransDataIrreg.WL{i})
        ki = dispersion_free_surface((2*pi*TransDataIrreg.F{i}(j))^2/Param.g,0,H);
        TransDataIrreg.WLN{i}(j,1) = 2*pi/imag(ki);
    end
    
end



figure()
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};
plot(WL,MultTDS,'-r')
PlotNames{end +1} = 'Wiki Waves';
hold on;
plot(WL,MultTLB,'-k')
PlotNames{end +1} = 'Lukes';
% plot(WL,MultTLBNS,'-b')
% PlotNames{end +1} = 'Lukes No Surge';
% plot(WL,MultTLBLFNS,'-g')
% PlotNames{end +1} = 'Lukes Long Floe No Surge';
plot(WL,Boltz,'--k')
PlotNames{end +1} = 'Boltzman';
plot(WLReg,TransDataReg.trans,'.k','MarkerSize',12)
PlotNames{end +1} = 'Lukes Trans';
%Irregular
for j = 1: length(TransDataIrreg.TpTarg)
    plot(TransDataIrreg.WL{j}, TransDataIrreg.Trans{j} , Symbols{j},'Color',Cols{j} ,'MarkerSize', 16);
    pname = ['Irregular Conc  |  Tranmission Probes | Hs ', num2str(TransDataIrreg.HsTarg{j}),' [m] | Tp ', num2str(TransDataIrreg.TpTarg{j}),' [s]'];
    PlotNames{end +1} = pname;  
end
% plot(WLRegNew,TransDataRegNew.TransConcCal,'.k','MarkerSize',12)
% PlotNames{end +1} = 'New Trans';

% plot(WLRegNew,TransDataRegNew.TransRawConc,'.r','MarkerSize',12)
% PlotNames{end +1} = 'Raw Reflection/ Transmission';
% plot(WLRegNew,TransDataRegNew.TransCorConc,'^b','MarkerSize',12)
% PlotNames{end +1} = 'Corrected Reflection/ Transmission';

xlabel('Wavelength (m)')
ylabel('Transmission Coefficient (Energy)')
title('Multiple Transmission')
legend(PlotNames);









