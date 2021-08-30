%%% Code written by David Skene during his PhD.
%%% It models the water going onto the edges of the floating elastic plate.
%%% It then saves the data.

clear all
close all
clc

%currently configured to run multiple tests via these opening for loops.

% WavePers = 0.6:0.05:2;
% Amps = 0:0.01:0.08;

% WavePers = 0.95;
% Amps = 0:0.005:0.04;

WavePers = 0.95;
Amps = 0.03;
% %Amps = 0:0.005:0.04;
% 
% %WavePers = 0.6:0.05:2;
% % Amps = 0.03;
% %Amps = 0:0.005:0.04;

%select wave period (s)
for waveperiod=WavePers%[0.8,0.9,1.0]
    %select incident wave steepness
    th_res=100;
    SURGE=0;
    terms_grn=100;
    extra_pts=[];     
    rigid = 4; 
    n=100;
    N = 50;
    if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid); 
        Param = ModParam_def(Param,1,n,extra_pts,terms_grn,th_res); end

    freq=1/waveperiod; %wave frequency Hz
    omega=2*pi*freq; %angular frequency rad/s
    alpha=omega^2/Param.g; %scaled frequency
    k= dispersion_free_surface(alpha,0,Param.bed); %wavenumber
    g = Param.g;
    kr = imag(k);
    L = Param.floe_diam;

    xbeam=-L/2:L/((2*N + 1)-1):L/2;
    
    % out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',SURGE,0,1,1,1,'r',xbea
    out = fn_ElasticRaft2d('freq',1/waveperiod,Param,'disk',SURGE,0,1,1,0,'r',xbeam);
    T = (out(1).value);
    R = (out(2).value);
    EtaOut = out(3).value;
    uSG = out(4).value;
%     xxOut = out(4).value;
   

    figure();
    draughtratio=Param.rho/Param.rho_0; %archemedies submergence
    plateabove=Param.thickness*(1-draughtratio); %portion of plate above water
    platebelow=Param.thickness*draughtratio; %portion of plate below water

    for amp=Amps
        
        stramp = num2str(round(amp,6));
        strwp = num2str(round(waveperiod,6));
        FileN = strcat('/Volumes/Storage/PRSA-Analysis/OW-Matrix/test_',stramp,'_',strwp,'_nohits_fixedmotion.mat');

        if ~isfile(FileN)
            %get wavenumber of incident waves
            I=amp; %incident wave height m
%             xbeam = -L/2:L/N:L/2; %positon of beam
            factor=2; %just something to define the length of the beam. Note that beam has length 2L.

            x_left=-factor*L/2:L/2/N:-L/2; %points to left of plate
            x_right=L/2:L/2/N:factor*L/2; %points to right of plate

            xp = @(x,t) (x -real(I*uSG*exp(-1i*omega*t)));
            leftwave=@(x,t) I*real( (exp(k*(xp(x,t) +L/2)) + R*exp(-k*(xp(x,t)+L/2)))*exp(-1i*omega*t));
            incwave=@(x,t) I*real( exp(k*((x)+L/2))*exp(-1i*omega*t));
            rightwave=@(x,t)I*real(T*exp(k*(xp(x,t)-L/2))*exp(-1i*omega*t));
            beamheight=@(j,t) I*real( EtaOut(j)*exp(-1i*omega*t))+plateabove;

            leftvelocity=@(x,t) I*real(g*k/(1i*omega)*( (exp(k*(xp(x,t)+L/2)) + R*exp(-k*(xp(x,t)+L/2)))*exp(-1i*omega*t)  )) + I*g*k/(1i*omega)*real(uSG*exp(-1i*omega*t)) ;
            rightvelocity=@(x,t) I*real(g*k/(1i*omega)*(T*exp(k*(xp(x,t)-L/2))*exp(-1i*omega*t))) + I*g*k/(1i*omega)*real(uSG*exp(-1i*omega*t));
            incvelocity=@(x,t) I*real(g*k/(1i*omega)*( (exp(k*(x+L/2)))*exp(-1i*omega*t) ));
    %           xxPL = -2.5:0.01:2.5
    %           figure()
    % 
    %           h1 = subplot(2,1,1); hold on; set(h1,'box','on')
    %           h2 = subplot(2,1,2); hold on; set(h2,'box','on')
    % 
    %           Incident = exp(k*(xxPL+L/2));
    % 
    %           IRreal = leftwave(xxPL(xxPL<-L/2),0);
    %           TIreal = rightwave(xxPL(xxPL>L/2),0);
    % 
    %           IRimag = leftwave(xxPL(xxPL<-L/2),pi/(2*omega));
    %           TIimag = rightwave(xxPL(xxPL>L/2),pi/(2*omega));
    %           
    %           plot(h1,xxPL,real(Incident),'k:')
    %           plot(h1,xxPL(xxPL<-L/2),IRreal,'b-')
    %           plot(h1,xxPL(xxPL>L/2),TIreal,'g-')
    %           plot(h1,xxOut,real(EtaOut),'-r')
    %           ylabel(h1,'Re(\eta)','fontsize',14)
    %           title(h1,['floe profile (r) and incident wave (k:), left wave (b-), right wave(g-)'],'fontsize',14)
    %           plot(h2,xxPL,imag(Incident),'k:')
    %           plot(h2,xxPL(xxPL<-L/2),IRimag,'b-')
    %           plot(h2,xxPL(xxPL>L/2),TIimag,'g-')
    %           plot(h2,xxOut,imag(EtaOut),'-r')
    %           xlabel(h2,'x','fontsize',14); ylabel(h2,'Im(\eta)','fontsize',14)

            dx=0.0005; %mesh size for SWEs

            beamleft=min(xbeam); %beam leftmost pos
            beamright=max(xbeam); %beam rightmost pos

            x=[beamleft:dx:beamright];

            %defining an initial time step. dx/5 seemed to keep things stable.
            dt=dx/5

            time=0;

            %defining a very small amount of water depth on the plate (the current SWE
            %solver does not do the wetting problem).
            fakemin=0.00001;

            %vector for the SWE variables. u(1,:) is the depth of the water. u(1,:) is
            %the depth multipled by the velocity.
            u=zeros(2,length(x));
            u(1,:)=u(1,:)+fakemin;

            %function handle for the SWE boundary conditions. Note that htis work
            %assumes a zero velocity at the plate edges for SWE forcing.
            uleft=@(t) [leftwave(max(beamleft),t);0];
            uright=@(t) [rightwave(min(beamright),t);0];

            %term for spatial discretisation
            spaceterm=u*0;

            %%%various terms for plotting the motion of the plate
            plotinterval=waveperiod/30;
            interval=plotinterval;
            % interval = 0;

            index=@(X) N/L*X+1+N;                       
            dipheight=@(t,index) beamheight(floor(index),t);

            %defining some initial conditions
%             xbeam=-L/2:L/(length(EtaOut)-1):L/2;
            u(:,1)=uleft(time)-[beamheight(1,time);0];
            u(:,length(u(1,:)))=uright(time)-[beamheight(length(EtaOut),time);0];

            u(1,1)=max(u(1,1),fakemin);
            u(1,length(u(1,:)))=max(u(1,length(u(1,:))),fakemin);

            %defining things for recording the solution.
            recordinterval=waveperiod/125;
            recordstartime=waveperiod*10;
            finaltime=5/waveperiod+recordstartime;
            rinterval=0;
            Hrecord=zeros(5000,length(u(1,:)));
            UHrecord=zeros(5000,length(u(1,:)));
            Trecord=zeros(5000,1);
            Beamrecord=zeros(5000,length(EtaOut));

            Vrecord=zeros(5000,1);
            VRightrecord=zeros(5000,1);
            Etaleftrecord=zeros(5000,1);
            Etarightrecord=zeros(5000,1);

            recordindex=0

            zlast=beamheight(1,time);
            zlastlast=beamheight(1,time);


            while time<finaltime

                %numerical solution to SWEs
                %this does the spatial discretisation
            %     [spaceterm,amax]=KTSpace_HLLE_mex(u,dx,x,1,fakemin*2);
                [spaceterm,amax]=KTSpace_HLLE_mex(u,dx,x,1,fakemin*2);

                %defines the time step, also makes sure it isnt too small.
                dt=min(1/amax/5*dx,0.00005);
                if dt<0.00005
                    dt=0.00005;
                end

                %using a simple first order euler for time step.
                u=u+dt*spaceterm;

                time=time+dt;

                %obtaining boundary conditions for next step
                 u(:,1)=uleft(time)-[beamheight(1,time);0];
                 u(1,1)=max(u(1,1),fakemin);

                 if u(1,1)>=fakemin*5
                     u(2,1)=max(u(1,1)*leftvelocity(beamleft,time),0);
                 end

                 u(:,length(u(1,:)))=(uright(time)-[beamheight(length(EtaOut),time);0] )   ;

                 %numerically putting a sink along the plate. Used because we assume
                 %energy is lost in overwash.
                  u(:,round(length(u)/2.5))=[0;0];
                 u(:,round(length(u)*(1-1/2.5)))=[0;0];

                 u(1,length(u(1,:)))=max(u(1,length(u(1,:))),fakemin); %fakemin

                 if u(1,length(u(1,:)))>=fakemin*5
                     u(2,end)=min(u(2,end)*rightvelocity(beamright,time),0);
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

                     Beamrecord(recordindex,:)=beamheight(1:length(EtaOut),time);
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

    %         %      plotting stuff
                 if interval>plotinterval
                     
%                     TopToInc =  beamheight(1:length(EtaOut),time) - incwave(xbeam,time).';
%                     [mean(TopToInc),plateabove]
    
                    subplot(2,2,2)
                    plot(x,u(1,:))
                    axis([min(x),max(x),0,0.02])
                    title('SWWE Overwash')
                    xlabel('x(m)')
                    ylabel('z(m)')
                    interval=0;
                    pause(0.01)
% %                     finaltime-time
% %                     steepness
% %                     waveperiod
    
                    subplot(2,2,1)
                    hold off
                    plot(x_left,leftwave(x_left,time))
                    hold on
                    plot(x_right,rightwave(x_right,time))
                    xp = -L*4:0.1:L*4;
                    plot(xp,incwave(xp,time),'--k')
                    plot(xbeam,beamheight(1:length(EtaOut),time),'-')
                    plot(xbeam,beamheight(1:length(EtaOut),time) - Param.thickness,'-')
                    title('Linear Model')
                    xlabel('x(m)')
                    ylabel('z(m)')
                    axis([min(x_left),max(x_right),-2*I- Param.thickness,2*I + Param.thickness])
                    
                    subplot(2,2,3)
                    hold off
                    plot(x_left,leftvelocity(x_left,time))
                    hold on
                    plot(x_right,rightvelocity(x_right,time))
                    plot(xp,incvelocity(xp,time),'--k')
                    title('Linear Model')
                    xlabel('x(m)')
                    ylabel('u(m/s)')
                    axis([min(x_left),max(x_right),-2*I*g*imag(k)/(omega)-10^-3,2*I*g*imag(k)/(omega)+10^-3])

                    subplot(2,2,4)
                    plot(x,u(2,:))
                    title('SWWE Overeash Model')
                    xlabel('x(m)')
                    ylabel('u(m/s)')
                    axis([min(x),max(x),-0.01,0.01])
                    interval=0;
                 end
             end
            %save result
            save(FileN)

            clearvars -except xbeam uSG WavePers Amps waveperiod th_res SURGE terms_grn extra_pts rigid n N Param freq omega alpha ki kr g k KAmp out T R EtaOut xxOut L draughtratio plateabove platebelow
        end
    end
end










