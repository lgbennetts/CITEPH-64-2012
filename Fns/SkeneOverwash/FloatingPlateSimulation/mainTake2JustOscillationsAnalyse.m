%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

% clear all
close all
% clc

%currently configured to run multiple tests via these opening for loops.

%select wave period (s)
for waveperiod=1/(2*pi)%0.95%[0.8,0.9,1.0]
%select incident wave steepness
for steepness= 0.14%4.4519

    clearvars -except steepness waveperiod

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

draughtratio=rho_i/rho_w ;%archemedies submergence
plateabove=thickness*(1-draughtratio); %portion of plate above water
platebelow=thickness*draughtratio; %portion of plate below water

%configuring some nondimensional variables relevant to LPT solution
gamma=rho_i*thickness/rho_w;
beta=E*thickness^3/(12*(1-nu^2)*rho_w*g);

D = E*thickness^2/(12*rho_i*(1-nu^2));
% CL = L;
% alpha = CL*omega^2/g;
% beta = (rho_i*thickness*D) / (rho_w*g*CL^4);
% gamma = rho_i*thickness/(rho_w*CL);

CL = L;
alpha = omega^2/g;
beta = (rho_i*thickness*D) / (rho_w*g);
gamma = rho_i*thickness/(rho_w);



1/sqrt(alpha);
g/omega/sqrt(alpha);
omega;

%Get linear potential theory solution
%displacement is plate's frequency domain displacement
%R and T are the reflected and transmitted complex amplitudes that are
%scaled
[xi,displacement,potential,R,T,f] = elastic_plate_modes(alpha,beta,gamma,H,L/2,nDS,N);
xbeam=-L/2:L/(length(displacement)-1):L/2;

[abs(R), abs(T)]
out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',SURGE,0,1,1,0,'r',xbeam);
TN = (out(1).value);
RN = (out(2).value);
EtaOut = out(3).value;
uSG = out(4).value;
% xxOut = out(4).value;
[abs(RN),abs(R), abs(TN),abs(T)]

%get wavenumber of incident waves
k=k(1);

I=steepness/(-1i*k); %incident wave height m

%no idea what this is. I believe I used it in debugging
T1k=abs(T*k)*I;


factor=4; %just something to define the length of the beam. Note that beam has length 2L.

x_left=-factor*L/2:L/2/N:-L/2; %points to left of plate
x_right=L/2:L/2/N:factor*L/2; %points to right of plate

%Using David's
%function handle for waves to the left and right of the plate.
leftwave=@(x,t) I*real(exp(-k*(-x))*exp(1i*omega*(-t)) + R*exp(k*(-x))*exp(1i*omega*(-t)));
rightwave=@(x,t) I*real(T*exp(-k*(-x))*exp(1i*omega*(-t)));
incwave=@(x,t) I*real(exp(-k*(-x))*exp(1i*omega*(-t)));
%function handle for the displacement of the ith point along the beam.
% beamheight=@(i,t) I*g/omega*real(1i/sqrt(alpha)*displacement(i)*exp(1i*omega*t))+plateabove;
% beamheight=@(i,t) I*real(1i/sqrt(alpha)*displacement(i)*exp(1i*omega*t))+plateabove;
beamheight=@(i,t) I*real(1i/sqrt(alpha)*displacement(end + 1 - i)*exp(1i*omega*(-t)))+plateabove;

%function handle for the velocity of the water at the left and right of the
%plate.
leftvelocity=@(x,t) I*real(-1i*g*k/omega*(exp(-k*(-x))*exp(1i*omega*t)+R*exp(k*(-x))*exp(1i*omega*(-t))));
rightvelocity=@(x,t) I*real(-1i*g*k/omega*(T*exp(-k*(-x))*exp(1i*omega*(-t))));
incvelocity=@(x,t) I*real(-1i*g*k/omega*(exp(-k*(-x))*exp(1i*omega*(-t))));

dx = xbeam(2) - xbeam(1);
PlateBC = PlateBoundary(1i*displacement,potential,alpha,gamma,beta,dx);
% a = 1.2;
% xC = -10:20/10000:10;
% dx = xC(2) - xC(1);
%Tests - Free Surface
% y = sin(a*xC);
% y4a = a^4*sin(a*xC);
% y4n = FourthDeriv(y,dx);
% 
% figure();
% plot(xC,y4a);
% hold on;
% plot(xC,y4n);
% 
% 
% D4 = FourthDeriv(real(- 1i*sqrt(alpha)*displacement),dx);
% 
% figure();
% plot(xbeam,real(- 1i*sqrt(alpha)*displacement));
% hold on;
% plot(xbeam,D4);
% % PlateBC = PlateBoundary( real(- 1i*sqrt(alpha)*potential),real(potential),alpha,gamma,beta,dx);


% %Using Mine
% % xp = @(x,t) x - real(I*uSG*exp(-1i*omega*t));
% % leftwaveN=@(x,t) I*real( (exp(k*(xp(x,t) +L/2)) + RN*exp(-k*(xp(x,t)+L/2)))*exp(-1i*omega*t));
% % incwaveN=@(x,t) I*real( exp(k*(xp(x,t)+L/2))*exp(-1i*omega*t));
% % rightwaveN=@(x,t)I*real(TN*exp(k*(xp(x,t)-L/2))*exp(-1i*omega*t));
% % beamheightN=@(j,t) I*real( EtaOut(j)*exp(-1i*omega*t))+plateabove;
% % 
% % leftvelocityN=@(x,t) I*real(g*k/(1i*omega)*( (exp(k*(xp(x,t)+L/2)) + RN*exp(-k*(xp(x,t)+L/2)))*exp(-1i*omega*t)  )) ;
% % rightvelocityN=@(x,t) I*real(g*k/(1i*omega)*(TN*exp(k*(xp(x,t)-L/2))*exp(-1i*omega*t)));
% % incvelocityN=@(x,t) I*real(g*k/(1i*omega)*( (exp(k*(xp(x,t)+L/2)))*exp(-1i*omega*t) ));
% 
% leftwaveN=@(x,t) I*real( (exp(k*(x +L/2)) + RN*exp(-k*(x+L/2)))*exp(-1i*omega*t));
% incwaveN=@(x,t) I*real( exp(k*(x+L/2))*exp(-1i*omega*t));
% rightwaveN=@(x,t)I*real(TN*exp(k*(x-L/2))*exp(-1i*omega*t));
% beamheightN=@(j,t) I*real( EtaOut(j)*exp(-1i*omega*t))+plateabove;
% 
% leftvelocityN=@(x,t) I*real(g*k/(1i*omega)*( (exp(k*(x+L/2)) + RN*exp(-k*(x+L/2)))*exp(-1i*omega*t)  )) ;
% rightvelocityN=@(x,t) I*real(g*k/(1i*omega)*(TN*exp(k*(x-L/2))*exp(-1i*omega*t)));
% incvelocityN=@(x,t) I*real(g*k/(1i*omega)*( (exp(k*(x+L/2)))*exp(-1i*omega*t) ));

% dt = 0.01;
% xwhole = -factor*L/2: 2*factor*L/100:factor*L/2;
% time = 0;
% finaltime = 5*waveperiod;
% figure();
% set(gca,'FontSize',24) 
% while time<finaltime
% % hold off;
% % plot(x_left,leftwave(x_left,time),'-k','DisplayName','Ocean Wave - Wiki Waves')
% % hold on;
% % plot(xwhole,incwave(xwhole,time),'.k','DisplayName',' Incident Ocean Wave - Wiki Waves')
% % plot(x_right,rightwave(x_right,time),'-k','DisplayName','Ocean Wave - Wiki Waves')
% % plot(xbeam,beamheight(1:length(displacement),time),'-b','DisplayName','Top and bottom of floe - Wiki Waves')
% % plot(xbeam,beamheight(1:length(displacement),time)-thickness,'-b','DisplayName','Top and bottom of floe - Wiki Waves')
% % 
% % plot(x_left,leftwaveN(x_left,time),'-r','DisplayName','Wave - code') 
% % plot(x_right,rightwaveN(x_right,time),'-r','DisplayName','Wave - code')
% % plot(xwhole,incwaveN(xwhole,time),'.r','DisplayName',' Incident Ocean Wave - code')
% % plot(xbeam,beamheightN(1:length(displacement),time),'--g','DisplayName','Top and bottom of floe - code')
% % plot(xbeam,beamheightN(1:length(displacement),time) - thickness,'--g','DisplayName','Top and bottom of floe - code')
% % axis([min(x_left),max(x_right),-2*I,2*I])
% % title('Plate Motion')
% % xlabel('x(m)')
% % ylabel('z(m)')
% % legend()
% % pause(0.01)
% 
% % TopToInc =  beamheightN(1:length(EtaOut),time) - incwaveN(xbeam,time).';
% % [mean(TopToInc),plateabove]
% 
% 
% subplot(2,1,1)
% hold off
% plot(x_left,leftwave(x_left,time),'-b','LineWidth',2)
% hold on
% plot(xwhole,incwave(xwhole,time),':b','LineWidth',2)
% plot(xbeam ,beamheight(1:length(displacement),time),'-r','LineWidth',2)
% plot(x_right,rightwave(x_right,time),'-b','LineWidth',2)
% plot(xbeam ,beamheight(1:length(displacement),time)-thickness,'-r','LineWidth',2)
% axis([min(x_left),max(x_right),-2*I-thickness,2*I+thickness])
% xlabel('x(m)')
% ylabel('z(m)')
% title('Wiki Waves')
% % legend('Modified WaveField','Incident Wavefield','Floe')
% 
% subplot(2,1,2)
% hold off
% plot(x_left,leftwaveN(x_left,time),'-b','LineWidth',2)
% hold on
% plot(xwhole,incwaveN(xwhole,time),':b','LineWidth',2)
% plot(xbeam + real(I*uSG*exp(-1i*omega*time)),beamheightN(1:length(displacement),time),'-r','LineWidth',2)
% plot(x_right,rightwaveN(x_right,time),'-b','LineWidth',2)
% plot(xbeam + real(I*uSG*exp(-1i*omega*time)),beamheightN(1:length(displacement),time)-thickness,'-r','LineWidth',2)
% axis([min(x_left),max(x_right),-2*I-thickness,2*I+thickness])
% xlabel('x(m)')
% ylabel('z(m)')
% title('Lukes')
% % legend('Modified WaveField','Incident Wavefield','Floe')
% 
% % subplot(2,2,3)
% % hold off
% % plot(x_left,leftvelocity(x_left,time),'-g','LineWidth',2)
% % hold on
% % plot(xwhole,incvelocity(xwhole,time),'--k','LineWidth',2)
% % plot(x_right,rightvelocity(x_right,time),'-g','LineWidth',2)
% % axis([min(x_left),max(x_right),min(-2*I*g*k/(1i*omega),-0.001),max(2*I*g*k/(1i*omega),0.001)])
% % title('Wiki Waves')
% % xlabel('x(m)')
% % ylabel('u(m/s)')
% % % legend('Modified WaveField','Incident Wavefield')
% % 
% % subplot(2,2,4)
% % hold off
% % plot(x_left,leftvelocityN(x_left,time),'-g','LineWidth',2)
% % hold on
% % plot(xwhole,incvelocityN(xwhole,time),'--k','LineWidth',2)
% % plot(x_right,rightvelocityN(x_right,time),'-g','LineWidth',2)
% % axis([min(x_left),max(x_right),min(-2*I*g*k/(1i*omega),-0.001),max(2*I*g*k/(1i*omega),0.001)])
% % xlabel('x(m)')
% % ylabel('u(m/s)')
% % title('Lukes')
% % legend('Modified WaveField','Incident Wavefield')
% sgtitle(['Comparison Tp = ',num2str(waveperiod), ' Wavelength = ',num2str(2*pi/imag(k)), ' Amplitude = ',num2str(I), ' Steepness = ',num2str(steepness),' t = ',num2str(time)])
% pause(0.01)
% 
% time = time + dt;
% 
% end



end
end

function D4 = FourthDeriv(fx,dx)
D4 = zeros(size(fx));
for i = 1:2
D4(i) = (3*fx(i)-14*fx(i+1)+26*fx(i+2)-24*fx(i+3)+11*fx(i+4) -2*fx(i+5)) / dx^4;
D4(end + 1 -i) = (3*fx(end + 1 -i)-14*fx(end + 1 -(i+1))+26*fx(end + 1 -(i+2))-24*fx(end + 1 - (i+3))+11*fx(end + 1- (i+4)) - 2*fx(end + 1- (i+5))) / dx^4;
end
%
for i = 3:length(fx)-2
    D4(i) = (1*fx(i-2)-4*fx(i-1)+6*fx(i)-4*fx(i+1)+1*fx(i+2)) / dx^4;
end
end


function LHS =  PlateBoundary(DzPotential,Potential,alpha,gamma,beta,dx)
 
DzTerms = (beta*FourthDeriv(DzPotential,dx)) - (gamma*alpha-1)*DzPotential;
Rest = -alpha*Potential;

LHS =  DzTerms + Rest;
end








