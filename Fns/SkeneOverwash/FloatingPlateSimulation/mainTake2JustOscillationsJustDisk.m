%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

% clear all
close all
% clc

%currently configured to run multiple tests via these opening for loops.

%select wave period (s)
for waveperiod=10%[0.8,0.9,1.0]
%select incident wave steepness
for steepness= 0.14%4.4519

    clearvars -except steepness waveperiod

 th_res=100;
SURGE=0;
terms_grn=100;
extra_pts=[];     
rigid = 4; 
N=50; %determines the number of points along beam (i.e. floating plate)
n=1; %number of eigenvalues
if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid); 
    Param = ModParam_def(Param,1,n,extra_pts,terms_grn,th_res); end
%  
    
E=Param.E; %plate elasticity Pa
rho_i=Param.rho; %plate density kg/m^3
nu=Param.nu; %poisson's ratio

thickness=Param.thickness; %plate thinkcness m

rho_w=Param.rho_0; %water density kg/m^3
freq=1/waveperiod; %wave frequency Hz
g=9.81; %gravitational accel m/s^2
omega=2*pi*freq; %angular frequency rad/s
H=Param.bed; %water depth m
L=Param.floe_diam; %beam length m
% N=800 %determines the number of points along beam (i.e. floating plate)
% n=400 %number of eigenvalues

alpha=omega^2/g; %scaled frequency
k=dispersion_free_surface(alpha,0,H); %wavenumber
kr = imag(k);

draughtratio=rho_i/rho_w ;%archemedies submergence
plateabove=thickness*(1-draughtratio); %portion of plate above water
platebelow=thickness*draughtratio; %portion of plate below water

%configuring some nondimensional variables relevant to LPT solution
gamma=rho_i*thickness/rho_w;
beta=E*thickness^3/(12*(1-nu^2)*rho_w/g);

1/sqrt(alpha);
g/omega/sqrt(alpha);
omega;

%Get linear potential theory solution
%displacement is plate's frequency domain displacement
%R and T are the reflected and transmitted complex amplitudes that are
%scaled
[xi,displacement,potential,R,T,f] = elastic_plate_modes(alpha,beta,gamma,H,L/2,n,N);
xbeam=-L/2:L/(length(displacement)-1):L/2;

out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',SURGE,0,1,1,1,'r',xbeam);
TN = (out(1).value);
RN = (out(2).value);
EtaOut = out(3).value.';
uSG = out(4).value;
Roots = out(5).value;
% xxOut = out(4).value;

%get wavenumber of incident waves
k=k(1);

I=steepness/(-1i*k); %incident wave height m

%no idea what this is. I believe I used it in debugging
T1k=abs(T*k)*I;

factor=2; %just something to define the length of the beam. Note that beam has length 2L.

x_left=-factor*L/2:L/2/N:-L/2; %points to left of plate
x_right=L/2:L/2/N:factor*L/2; %points to right of plate

%Using David's
%function handle for waves to the left and right of the plate.
incwave=@(x,t) I*real(exp(-k*(x))*exp(1i*omega*t));
%function handle for the displacement of the ith point along the beam.
beamheight=@(i,t) I*g/omega*real(1i/sqrt(alpha)*displacement(i)*exp(1i*omega*t)) + plateabove;


%Using Mine
% ScaledEtaOut = EtaOut(end:-1:1) ./ abs(EtaOut(end:-1:1));

% beamheightN=@(j,t) I*real( EtaOut(end + 1 - j)*exp(1i*omega*(-t + kr/omega*(L/2) )))+plateabove;


beamheightN=@(j,t) I*real( (EtaOut(j)) *exp(-1i*omega*t))+ plateabove;
incwaveN=@(x,t) I*real(exp(k*(x+ L/2))*exp(-1i*omega*t));


dt = 0.01;
xwhole = -2*pi/(2*kr):0.01 :2*pi/(2*kr);
time = 0;
finaltime = 5*waveperiod;
figure();
set(gca,'FontSize',24) 
while time<finaltime
% hold off;
% plot(x_left,leftwave(x_left,time),'-k','DisplayName','Ocean Wave - Wiki Waves')
% hold on;
% plot(xwhole,incwave(xwhole,time),'.k','DisplayName',' Incident Ocean Wave - Wiki Waves')
% plot(x_right,rightwave(x_right,time),'-k','DisplayName','Ocean Wave - Wiki Waves')
% plot(xbeam,beamheight(1:length(displacement),time),'-b','DisplayName','Top and bottom of floe - Wiki Waves')
% plot(xbeam,beamheight(1:length(displacement),time)-thickness,'-b','DisplayName','Top and bottom of floe - Wiki Waves')
% 
% plot(x_left,leftwaveN(x_left,time),'-r','DisplayName','Wave - code') 
% plot(x_right,rightwaveN(x_right,time),'-r','DisplayName','Wave - code')
% plot(xwhole,incwaveN(xwhole,time),'.r','DisplayName',' Incident Ocean Wave - code')
% plot(xbeam,beamheightN(1:length(displacement),time),'--g','DisplayName','Top and bottom of floe - code')
% plot(xbeam,beamheightN(1:length(displacement),time) - thickness,'--g','DisplayName','Top and bottom of floe - code')
% axis([min(x_left),max(x_right),-2*I,2*I])
% title('Plate Motion')
% xlabel('x(m)')
% ylabel('z(m)')
% legend()
% pause(0.01)


hold off
plot(xbeam,beamheight(1:length(displacement),time),'-r','LineWidth',2)
hold on
plot(xbeam,beamheight(1:length(displacement),time)-thickness,'-r','LineWidth',2)
plot(xbeam,beamheightN(1:length(displacement),time),'-b','LineWidth',2)
plot(xbeam,beamheightN(1:length(displacement),time)-thickness,'-b','LineWidth',2)
% plot(xwhole, incwave(xwhole,time),'-k','LineWidth',2)
plot(xwhole, incwaveN(xwhole,time),'--k','LineWidth',2)
axis([min(xwhole),max(xwhole),-2*I-thickness,2*I+thickness])
xlabel('x(m)')
ylabel('z(m)')
title(['Comparison - Tp = ',num2str(waveperiod), ' WaveLength = ',num2str(2*pi/kr), 't = ',num2str(time)])
pause(0.01)


time = time + dt;

end



end
end










