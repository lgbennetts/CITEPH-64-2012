%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

% clear all
close all
% clc

SymsCols1= {'-b','-r','-k','-g','--b','--r'};
SymsCols2= {'-.b','-.r','-.k','-.g',':b',':r'}; 

NVals = [200];
WavePer = 0.5:0.1:2;

TDS = zeros(length(NVals),length(WavePer));
TLB = zeros(length(NVals),length(WavePer));
WL = zeros(length(NVals),length(WavePer));

th_res=100;
extra_pts=[];
terms_grn=100;


for j = 1:length(NVals)
%currently configured to run multiple tests via these opening for loops.
SURGE=0;
rigid = 100; 
N=NVals(j); %determines the number of points along beam (i.e. floating plate)
n = NVals(j);
Param = ParamDef_Oceanide(rigid); 
Param = ModParam_def(Param,1,N,extra_pts,terms_grn,th_res);
%  
    
E=Param.E; %plate elasticity Pa
rho_i=Param.rho; %plate density kg/m^3
nu=Param.nu; %poisson's ratio

thickness=Param.thickness; %plate thinkcness m
rho_w=Param.rho_0; %water density kg/m^3
g=Param.g; %gravitational accel m/s^2
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

    D = E*thickness^2/(12*rho_i*(1-nu^2));
    CL = L;
    alpha = CL*omega^2/g;
    beta = (rho_i*thickness*D) / (rho_w*g*CL^4);
    gamma = rho_i*thickness/(rho_w*CL);
    
    %Get linear potential theory solution
    %displacement is plate's frequency domain displacement
    %R and T are the reflected and transmitted complex amplitudes that are
    %scaled
    [xi,displacement,potential,R,T,f] = elastic_plate_modes(alpha,beta,gamma,H,L/2,n,N);
    xbeam=-L/2:L/(length(displacement)-1):L/2;

    out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',1,0,0,1,0,'r',xbeam);
    TN = (out(1).value);
    RN = (out(2).value);

    TDS(j,i) = abs(T)^2;
    TLB(j,i) = abs(TN)^2;
    WL(j,i) = 2*pi/kr;
    


end
end

%basic
figure()
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};
for j = 1:length(NVals)
    plot(WL(j,:),TDS(j,:),SymsCols1{j},'LineWidth',2)
    PlotNames{end +1} = ['Wiki Waves ', num2str(NVals(j))];
    hold on;
    plot(WL(j,:),TLB(j,:),SymsCols2{j},'LineWidth',2)
    PlotNames{end +1} = ['Lukes ', num2str(NVals(j))];
end
xlabel('Wavelength (m)')
ylabel('Transmission Coefficient (Energy)')
title('Convergence')
legend(PlotNames);






