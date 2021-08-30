%%% Code written by David Skene during his PhD.
%%% It loads the prior simulations and calculates the corrected wave
%%% transmission via the transition loss theory method of Skene & Bennetts 2021

clear all
close all
clc

%just opening all the experiment plots at once. Did this because I was a
%PhD student.
% uiopen('.\T08tsm_noblack.fig',1)
% uiopen('.\T09tsm_noblack.fig',1)
% uiopen('.\T10tsm_noblack.fig',1)


WavePers = 0.6 %0.3:0.05:0.6;


th_res=100;
SURGE=1;
terms_grn=100;
extra_pts=[];
rigid = 4; 
n=100;
N = 50;
if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid); 
    Param = ModParam_def(Param,1,n,extra_pts,terms_grn,th_res); end

%select period to calculate and plot results for
for periodi =1 : length(WavePers)
    
    periodname = num2str(round(WavePers(periodi),6));
    
    
    clearvars -except periodi Lfactor periodname Lrightfactor WavePers Param
    
    %vectors for old and new transmitted steepness and incident wave steepness.
    kT=[0]
    kTnew=[0]
    kI=[0]

    freq=1/WavePers(periodi); %wave frequency Hz
    omega=2*pi*freq; %angular frequency rad/s
    alpha=omega^2/Param.g; %scaled frequency
    k= dispersion_free_surface(alpha,0,Param.bed); %wavenumber
    g = Param.g;
    kr = imag(k);
    
    KAmp = 0:0.05:1;
    Amps = KAmp/kr;
    
    %running each individual steepness test.
    for ampi=1 : length(Amps) 

        ampname = num2str(round(Amps(ampi),6));

        %load that particular period and steepness simulation
        load(strcat('/Volumes/Storage/PRSA-Analysis/OW-Matrix/test_',ampname,'_',periodname,'_nohits_fixedmotion.mat')) 

        %do a little transform on the key saved variables to make furture code
        %easier.
        Hrecord=Hrecord';
        UHrecord=UHrecord';
        Trecord=Trecord-Trecord(1)

        %define linear potential theory reflection and transmission amplitudes.
        R=abs(R)*I;
        T=abs(T)*I;
        % abs(T);
        % 
        % T = Trans*I;
        % R = (1 - Trans)*I;

        %define energy of incident wave
        Ein=1/4*g^2/omega*(I^2-R^2);

        %This is the part where the code calculates the energy going into the
        %overwash.
        period=waveperiod;
        sampleperiods=2; %periods to average over
        tstart=period*4; %time to start averaging
        tend=tstart+period*sampleperiods; %time to end averaging
        [dummy1,dummy2]=max(Trecord); %removing zero entries in time vector
        Trecord=Trecord(1:dummy2);

        %getting starting and ending time indicies
        indextstart=max(find(Trecord<tstart));
        indextend=max(find(Trecord<tend));

        startpos=1 %first domain location
        srecord=[]
        probeA=3 %probe at leading edge of the plate, location A

        endpos=length(Hrecord(:,1))
        endprobe=endpos-3
        probeB=endprobe %probe at trailing edge of the plate, location B

        %Note: both of these are the third probe in. This is because location 1 is
        %the boundary condition. Location 2 is where the shock always is. And
        %therefore location 3 is where the upwave/downwave edge is actually
        %modelled in the solver.

        htest=[]
        Eplusavg=0

        %doing the process of calculating the energy going into the overwash at
        %each timestep.
        AvgEow=0
        AvgEAedge=0
        AvgEBedge=0
        for i=indextstart:indextend-1

            dt=Trecord(i+1)-Trecord(i);

            Azbeam=Beamrecord(i,1);
            Ahminus=Etaleftrecord(i);
            Ahplus=Hrecord(probeA,i);
            Auminus=Vrecord(i);
            Auplus=UHrecord(probeA,i)/Ahplus;
            if Ahplus<=fakemin
                Ahplus=0;
                Auminus=0;
            end
           if Ahminus<Azbeam
                Ahminus=0;
                Auminus=0;
           end


            Bzbeam=Beamrecord(i,length(EtaOut));
        %     Bzbeam=Beamrecord(i,length(displacement));
            Bhminus=Etarightrecord(i);
            Bhplus=Hrecord(probeB,i);
            Buminus=VRightrecord(i);
            Buplus=UHrecord(probeB,i)/Bhplus;
            if Bhplus<=fakemin
                Bhplus=0;
                Buminus=0;
            end
           if Bhminus<Bzbeam
                Bhminus=0;
                Buminus=0;
           end

           EAow = (1/2*Auplus^3*Ahplus+g*Ahplus^2*Auplus);
           EBow =-(1/2*Buplus^3*Bhplus+g*Bhplus^2*Buplus);

           EAedge = -g*Auminus*(Ahminus-Azbeam)*(Ahminus-Azbeam)*0
           EBedge = +g*Buminus*(Bhminus-Bzbeam)*(Bhminus-Bzbeam)*0

           AvgEAedge=AvgEAedge+(EAedge)*dt
           AvgEBedge=AvgEBedge+(EBedge)*dt
           AvgEow=AvgEow+2*EAow*dt+2*EBow*dt
        end
        %average the energy going over the selected sampling duration
        AvgEow=AvgEow/(Trecord(indextend)-Trecord(indextstart))
        AvgEAedge=AvgEAedge/(Trecord(indextend)-Trecord(indextstart))
        AvgEBedge=AvgEAedge/(Trecord(indextend)-Trecord(indextstart))

        omega=omega

        %do the calculation for the transition loss corrected wave transmission
        NEWT=sqrt(T^2-(AvgEow+AvgEAedge+AvgEBedge)*4*omega/g^2)

        %do stuff to plot it all.
        k=abs(k)

        'ka_i'
        kI=[kI,k*I]

        'ka_t'
        kT=[kT,k*T]

        'ka_tnew'
        kTnew=[kTnew,NEWT*k]

        wavelength=2*pi/k
    end

figure()
hold on
if periodi<100000000
plot(kI,kTnew,'k-','LineWidth',2)
plot(kI,kT,'r--','LineWidth',2)
else
plot(kI,kTnew,'mo-')
plot(kI,kT,'c-')
end

axis([0.00,1,0.00,1])
title(strcat('Wavelength: ',num2str(wavelength),' m'))
xlabel('Incident wave steepness')
ylabel('Transmitted wave steepness')
pause(0.01)
%legend('show')
%legend('Linear theory','Overwash theory')

%Save Data
%Save Matrices
% MatFile_NM = strcat('./Data/OW/Tp',periodname,'SingleFloeOW.mat');
% save(MatFile_NM,'kI','kTnew','kT');
end