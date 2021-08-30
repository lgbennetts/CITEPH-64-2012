%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

clear all
close all
clc

%currently configured to run multiple tests via these opening for loops.

%select wave period (s)
for waveperiod=0.9%[0.8,0.9,1.0]
%select incident wave steepness
for steepness=0.04%[0.04,0.06,0.08,0.1,0.12,0.14,0.15]

    clearvars -except steepness waveperiod
    
E=500e6 %plate elasticity Pa
rho_i=570 %plate density kg/m^3
nu=0.4 %poisson's ratio

thickness=10*10^-3 %plate thinkcness m

rho_w=1000 %water density kg/m^3
freq=1/waveperiod %wave frequency Hz
g=9.81 %gravitational accel m/s^2
omega=2*pi*freq %angular frequency rad/s
H=0.9 %water depth m
L=1.0 %beam length m
N=800 %determines the number of points along beam (i.e. floating plate)
n=400 %number of eigenvalues

alpha=omega^2/g %scaled frequency
k=dispersion_free_surface(alpha,0,H) %wavenumber

draughtratio=rho_i/rho_w %archemedies submergence
plateabove=thickness*(1-draughtratio) %portion of plate above water
platebelow=thickness*draughtratio %portion of plate below water

%configuring some nondimensional variables relevant to LPT solution
gamma=rho_i*thickness/rho_w
beta=E*thickness^3/(12*(1-nu^2)*rho_w/g)

1/sqrt(alpha)
g/omega/sqrt(alpha)
omega

%Get linear potential theory solution
%displacement is plate's frequency domain displacement
%R and T are the reflected and transmitted complex amplitudes that are
%scaled
[xi,displacement,potential,R,T,f] = elastic_plate_modes(alpha,beta,gamma,H,L/2,n,N);

%get wavenumber of incident waves
k=k(1)

I=steepness/(-1i*k) %incident wave height m

%no idea what this is. I believe I used it in debugging
T1k=abs(T*k)*I

xbeam = -L/2:L/N:L/2 %positon of beam

factor=2 %just something to define the length of the beam. Note that beam has length 2L.

x_left=-factor*L/2:L/2/N:-L/2; %points to left of plate
x_right=L/2:L/2/N:factor*L/2; %points to right of plate

%function handle for waves to the left and right of the plate.
leftwave=@(x,t) I*real(exp(-k*(x))*exp(1i*omega*t) + R*exp(k*(x))*exp(1i*omega*t));
rightwave=@(x,t) I*real(T*exp(-k*(x))*exp(1i*omega*t));

%function handle for the displacement of the ith point along the beam.
beamheight=@(i,t) I*g/omega*real(1i/sqrt(alpha)*displacement(i)*exp(1i*omega*t))+plateabove;

%function handle for the velocity of the water at the left and right of the
%plate.
leftvelocity=@(x,t) I*real(-1i*g*k/omega*(exp(-k*(x))*exp(1i*omega*t)+R*exp(k*(x))*exp(1i*omega*t)));
rightvelocity=@(x,t) I*real(-1i*g*k/omega*(T*exp(-k*(x))*exp(1i*omega*t)));



dx=0.0005 %mesh size for SWEs

beamleft=min(xbeam) %beam leftmost pos
beamright=max(xbeam) %beam rightmost pos

x=[beamleft:dx:beamright]

%defining an initial time step. dx/5 seemed to keep things stable.
dt=dx/5

time=0

%defining a very small amount of water depth on the plate (the current SWE
%solver does not do the wetting problem).
fakemin=0.00001

%vector for the SWE variables. u(1,:) is the depth of the water. u(1,:) is
%the depth multipled by the velocity.
u=zeros(2,length(x))
u(1,:)=u(1,:)+fakemin

%function handle for the SWE boundary conditions. Note that htis work
%assumes a zero velocity at the plate edges for SWE forcing.
uleft=@(t) [leftwave(max(beamleft),t);0]
uright=@(t) [rightwave(min(beamright),t);0];

%term for spatial discretisation
spaceterm=u*0;

%%%various terms for plotting the motion of the plate
plotinterval=waveperiod/30;
interval=plotinterval;

index=@(X) N/L*X+1+N;                       
dipheight=@(t,index) beamheight(floor(index),t);

%defining some initial conditions
xbeam=-L/2:L/(length(displacement)-1):L/2
u(:,1)=uleft(time)-[beamheight(1,time);0];
u(:,length(u(1,:)))=uright(time)-[beamheight(length(displacement),time);0];

u(1,1)=max(u(1,1),fakemin);
u(1,length(u(1,:)))=max(u(1,length(u(1,:))),fakemin);

%defining things for recording the solution.
recordinterval=waveperiod/125
recordstartime=waveperiod*10
finaltime=5/waveperiod+recordstartime
rinterval=0
Hrecord=zeros(5000,length(u(1,:)));
UHrecord=zeros(5000,length(u(1,:)));
Trecord=zeros(5000,1)
Beamrecord=zeros(5000,length(displacement));

Vrecord=zeros(5000,1)
VRightrecord=zeros(5000,1)
Etaleftrecord=zeros(5000,1)
Etarightrecord=zeros(5000,1)

recordindex=0

zlast=beamheight(1,time)
zlastlast=beamheight(1,time)
     
while time<finaltime
    
    %numerical solution to SWEs
    %this does the spatial discretisation
%     [spaceterm,amax]=KTSpace_HLLE_mex(u,dx,x,1,fakemin*2);
    [spaceterm,amax]=KTSpace_HLLE_mex(u,dx,x,1,fakemin*2);
    
    %defines the time step, also makes sure it isnt too small.
    dt=min(1/amax/5*dx,0.00005);
    if dt<0.00005
        dt=0.00005
    end
    
    %using a simple first order euler for time step.
    u=u+dt*spaceterm;

    time=time+dt;

    %obtaining boundary conditions for next step
     u(:,1)=uleft(time)-[beamheight(1,time);0];
     u(1,1)=max(u(1,1),fakemin);
     
     if u(1,1)>=fakemin*5
         u(2,1)=u(1,1)*leftvelocity(beamleft,time);
     end
     
     u(:,length(u(1,:)))=(uright(time)-[beamheight(length(displacement),time);0] )   ;

     %numerically putting a sink along the plate. Used because we assume
     %energy is lost in overwash.
     u(:,round(length(u)/2.5))=[0;0];
     u(:,round(length(u)*(1-1/2.5)))=[0;0];
     
     u(1,length(u(1,:)))=max(u(1,length(u(1,:))),fakemin); %fakemin
     
     if u(1,length(u(1,:)))>=fakemin*5
         u(2,end)=u(2,end)*rightvelocity(beamright,time);
     end
     
     %updating the intervals for recording purposes
     interval=interval+dt;
     rinterval=rinterval+dt;
     
     %recording stuff
     if rinterval>recordinterval && time>recordstartime
         recordindex=recordindex+1;
         Hrecord(recordindex,:)=u(1,:);
         UHrecord(recordindex,:)=u(2,:);
         
         Vrecord(recordindex)=leftvelocity(max(beamleft),time);
         VRightrecord(recordindex)=rightvelocity(min(beamright),time);
         Etaleftrecord(recordindex)=leftwave(max(beamleft),time);
         Etarightrecord(recordindex)=rightwave(min(beamright),time);

         Beamrecord(recordindex,:)=beamheight(1:length(displacement),time);
         Trecord(recordindex)=time;
         rinterval=0;
         
         
         
         Vminus(recordindex)=leftvelocity(max(beamleft),time);
         bestprobe=3;
         if u(1,bestprobe)>fakemin
            Vplus(recordindex)=u(2,bestprobe)/u(1,bestprobe);
         else
            Vplus(recordindex)=0;
         end
         
         Etaminus(recordindex)=leftwave(max(beamleft),time);
         Zbeam(recordindex)=beamheight(1,time);
         Etaplus(recordindex)=u(1,3)+beamheight(1,time);
         
         
     end     
     
     %plotting stuff
     if interval>plotinterval
         
         subplot(4,1,1)
          plot(x,u(1,:))
          axis([min(x),max(x),0,0.02])
         interval=0
         pause(0.01)
        finaltime-time
        steepness
        waveperiod
        
        subplot(4,1,2)
        hold off
        plot(x_left,leftwave(x_left,time))
        hold on
        plot(x_right,rightwave(x_right,time))
        plot(xbeam,beamheight(1:length(displacement),time))
        axis([min(x_left),max(x_right),-2*I,2*I])
        
        subplot(4,1,3)
        hold off
        plot(x_left,leftvelocity(x_left,time))
        axis([min(x_left),max(x_right),-2*I*abs(k),2*I*abs(k)])
        hold on
        plot(x_right,rightvelocity(x_right,time))
        
%         subplot(4,1,4)
%         hold on
%             plot(time,bvz*omega/g,'ro')
%             plot(time,avz/g,'b*')
% %         plot(time,bvz*omega,'ro')
% %         plot(time,g,'k*')
     end
%      zlastlast=zlast
%      zlast=znow
     
end


%save result
% save(strcat('test_',num2str(steepness),'_',num2str(waveperiod),'_nohits_fixedmotion.mat'))


end
end










