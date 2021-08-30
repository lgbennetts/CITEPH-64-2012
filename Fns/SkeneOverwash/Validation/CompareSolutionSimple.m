%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

clear all
close all
% clc

SymsCols1= {'-b','-r','-k','-g','--b','--r'};
SymsCols2= {'-.b','-.r','-.k','-.g',':b',':r'}; 

NVals = 200;
WavePer = 0.95;%2*pi;%0.5:0.1:2;

TDS = zeros(length(NVals),length(WavePer));
TLB = zeros(length(NVals),length(WavePer));
WL = zeros(length(NVals),length(WavePer));

th_res=100;
extra_pts=[];
terms_grn=100;
rigid = 4;

I = 0.03;
%currently configured to run multiple tests via these opening for loops.
N=NVals; %determines the number of points along beam (i.e. floating plate)
n = NVals;
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
%select incident wave steepness
waveperiod=WavePer;

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

out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',0,0,0,1,0,'r',xbeam);
TN = (out(1).value);
RN = (out(2).value);
EtaOut = (out(3).value);

% TDS(j,i) = abs(T)^2;
% TLB(j,i) = abs(TN)^2;
% WL(j,i) = 2*pi/kr;
    

leftwave=@(x,t) I*real(exp(-k*x)*exp(1i*omega*t) + R*exp(k*(-x))*exp(1i*omega*t));
rightwave=@(x,t) I*real(T*exp(-k*x)*exp(1i*omega*t));
incwave=@(x,t) I*real(exp(-k*x)*exp(1i*omega*t));
beamheight=@(i,t) I*real(1i/sqrt(alpha)*displacement(i)*exp(1i*omega*t))+plateabove;

leftwaveN=@(x,t) I*real( (exp(k*(x +L/2)) + RN*exp(-k*(x+L/2)))*exp(-1i*omega*t));
incwaveN=@(x,t) I*real( exp(k*(x+L/2))*exp(-1i*omega*t));
rightwaveN=@(x,t)I*real(TN*exp(k*(x-L/2))*exp(-1i*omega*t));
beamheightN=@(j,t) I*real( (EtaOut(j))*exp(-1i*omega*t))+plateabove;

dt = WavePer/100.0;
factor=4; %just something to define the length of the beam. Note that beam has length 2L.
x_left=-factor*L/2:L/2/N:-L/2; %points to left of plate
x_right=L/2:L/2/N:factor*L/2; %points to right of plate
xwhole = -factor*L/2: 2*factor*L/100:factor*L/2;
time = 0;
finaltime = 10*waveperiod;
figure();
set(gca,'FontSize',24) 
while time<finaltime
    
    subplot(2,1,1)
    hold off
    plot(x_left,leftwave(x_left,time),'-b','LineWidth',2)
    hold on
    plot(xwhole,incwave(xwhole,time),':b','LineWidth',2)
    plot(xbeam ,beamheight(1:length(displacement),time),'-r','LineWidth',2)
    plot(x_right,rightwave(x_right,time),'-b','LineWidth',2)
    plot(xbeam ,beamheight(1:length(displacement),time)-thickness,'-r','LineWidth',2)
    axis([min(x_left),max(x_right),-2*I-thickness,2*I+thickness])
    xlabel('x(m)')
    ylabel('z(m)')
    title('Wiki Waves')
    % legend('Modified WaveField','Incident Wavefield','Floe')

    subplot(2,1,2)
    hold off
    plot(x_left,leftwaveN(x_left,time),'-b','LineWidth',2)
    hold on
    plot(xwhole,incwaveN(xwhole,time),':b','LineWidth',2)
    plot(xbeam,beamheightN(1:length(displacement),time),'-r','LineWidth',2)
    plot(x_right,rightwaveN(x_right,time),'-b','LineWidth',2)
    plot(xbeam,beamheightN(1:length(displacement),time)-thickness,'-r','LineWidth',2)
    axis([min(x_left),max(x_right),-2*I-thickness,2*I+thickness])
    xlabel('x(m)')
    ylabel('z(m)')
    title('Lukes')
    sgtitle(['Comparison Tp = ',num2str(waveperiod), ' Wavelength = ',num2str(2*pi/imag(k)), ' Amplitude = ',num2str(I),' t = ',num2str(time)])

    pause(0.01)

    
    time = time + dt;

end








