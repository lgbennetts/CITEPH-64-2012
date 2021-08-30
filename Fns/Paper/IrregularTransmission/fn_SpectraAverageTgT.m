% function fn_InvestigateSpectra
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
%
% General:
%
% conc      = either 39 or 79 , otherwise calibration
% WaveType =  'Regular' or 'Irregular'
% TestName = name of test in fn_*_testspecs function for asociated
% experiment
% probe - probe to read, allow any from  1- 20 (10 on left, and 10 on right of probe)
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function out = fn_SpectraAverageTgT(Tp,conc,TestName,WaveType,probes,OutputStr,wbin,TpersStr,t0description,t1description)

% 
if ~exist('conc','var') conc=79; end; %39,79 or empty (1)
if ~exist('probes','var'); probes= 11:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= '19'; end ;  
if ~exist('WaveType','var'); WaveType= 'Regular'; end ;  

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end;
if ~exist('wbin','var'); wbin =3; end ;  

if ~exist('t0description','var'); t0description ='waves reach x'; end ; 
% if ~exist('t1description','var'); t1description ='final waves reach x'; end ; 

if ~exist('t1description','var')
    if strcmp(WaveType,'Regular')
        if probes(1) == 1
            t1description ='MIZ ref waves reach x';
        else
            t1description ='beach ref waves reach x';
        end
        t1description ='final waves reach x';
    else
        t1description ='final waves reach x';
    end
end 



if ~exist('OutputStr','var'); OutputStr ='TargetSpectra Mean_AllProbe_Spectra Ind_Spectra_Whole'; end ;  

% close all;

out_str = ' ''dummy'' '; out_val = ' 0 ';

for i = 1:length(probes)
 %Read Data in Calib 1, at probe
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probes(i));
    
    % Tp = c_pram.period/10.0;
    if ~exist('Tp','var'); Tp = c_pram.period/10.0; end ;  
    
    dt = tm(2) - tm(1); 
    

    eval(TpersStr);
    T_windowSec = round(Tpers*Tp);
    SamplingRate = floor(1./dt);
    T_window = SamplingRate*T_windowSec;

    Hs = c_pram.wave_height/100.0;

    Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower

    t0index = find(strcmp({Tind.description},t0description));
    t0 = min(Tind(t0index).time,tm(end)); %+ T_windowSec/2;
    t1index = find(strcmp({Tind.description},t1description));
    t1 = min(Tind(t1index).time,tm(end));% - T_windowSec/2 ;
    
%     figure('DefaultAxesFontSize',18);
%      plot(tm,disp,'-k', 'DisplayName','Time Series Calibration Probe 11');
%      hold on;
%      plot([0,300],[Hs/2 Hs/2],'--r', 'DisplayName','Hs/2');
%      plot([0,300],[-Hs/2 -Hs/2],'--r', 'DisplayName','-Hs/2');
%       plot([t0,t0],[-0.025 0.025],'--b', 'DisplayName','t0');
%       plot([t1,t1],[-0.025 0.025],'--b', 'DisplayName','t1');
%      axis([0 300 -0.025 0.025 ])
%      title(['Elevation Time Series for Irregular Experiment with \tau = ', num2str(Tp), ' Hs = ', num2str(Hs)]);
%      xlabel('time (s)')
%      ylabel('elevation (m)')
%      legend();
        
    
    [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
    period = 1./ff;
    FourierWindow_Sec = FourierWindow/SamplingRate;

    jj0=find(tm_vec>t0,1);
    jj1=find(tm_vec<t1,1,'last');

    Sn1(:,i) = mean(Sn_mat(:,jj0:jj1), 2);
    
    
    if strfind(OutputStr,'Ind_Spectra_Whole')
        %Write a matrix for each
        SnA{i} = Sn_mat;
    end
    
end

MeanSpectra  = mean(Sn1,2);


%Outputs
if strfind(OutputStr,'Ind_Spectra_Whole')
 out_str = [out_str '; ''SpectralMat'' '];
 out_val = [out_val '; SnA'];
 out_str = [out_str '; ''PeriodMat'' '];
 out_val = [out_val '; period'];
end

if strfind(OutputStr,'Mean_AllProbe_Spectra')
 out_str = [out_str '; ''AllProbesMeanSpectra'' '];
 out_val = [out_val '; Sn1'];
 
 out_str = [out_str '; ''MeanSpectra'' '];
 out_val = [out_val '; MeanSpectra'];
 out_str = [out_str '; ''Period'' '];
 out_val = [out_val '; period'];
end


if strfind(OutputStr,'TargetSpectra')
 [~,jj]=min(abs(period-Tp));
  PerFilt = Tp;
  SpecFilt =(1./(period(jj - wbin) - period(jj + wbin))).* trapz(period(jj + [wbin:-1:-wbin]),MeanSpectra(jj + [wbin:-1:-wbin]));

%    figure();
%     plot(period, sqrt(2*MeanSpectra))
%     hold on;
%     plot([period(jj-wbin),period(jj-wbin)],[0,max(sqrt(2*MeanSpectra))],'--r', 'DisplayName','lim');
%     plot([period(jj+wbin),period(jj+wbin)],[0,max(sqrt(2*MeanSpectra))],'--r', 'DisplayName','lim');
%     plot([period(jj-wbin),period(jj+wbin)],[sqrt(2*SpecFilt),sqrt(2*SpecFilt)],'--k', 'DisplayName','average');
%    
%     xlabel('period (s)')
%      ylabel('a(f) (m)')

%sum(MeanSpectra(jj + [-wbin:wbin]));
%  SpecFilt = sum(MeanSpectra(jj + [-wbin:wbin]));
 
 out_str = [out_str '; ''PerFilt'' '];
 out_val = [out_val '; PerFilt'];
 out_str = [out_str '; ''SpecFilt'' '];
 out_val = [out_val '; SpecFilt'];
end

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});']) 
out(1)=[];

return

function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

