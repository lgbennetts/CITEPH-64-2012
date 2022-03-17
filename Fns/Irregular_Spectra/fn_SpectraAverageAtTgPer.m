% function fn_SpectraAverageAtTgPer
%
% DESCRIPTION: Extract spectra value at target period.
%         each spectra at target period is extracted from a single probe (average over frequency bins and between begin time and end time)
%         the spectract at target period is averaged over multiple probes 
%
% INPUTS:
%         Tp - Target period to extract spectra at
%         conc - Concentration of disk array
%         TestName - name of test (test_specifics)
%         WaveType - Regular or Irregular
%         probes - list of probes to extract from
%         OutputStr - String determining the outputs
%         wbin - size of frequency bin around Target to average over (length index)
%         TpersStr - String to evaluate the number of periods to extract
%         over (fn_Tpers)
%         t0description - string of description for start time 
%         t1description - string of description for end time
%
% OUTPUTS:
%     MatFile_NM - saved as a matrix file
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by

function out = fn_SpectraAverageAtTgPer(Tp,conc,TestName,WaveType,probes,OutputStr,wbin,TpersStr,t0description,t1description)

% 
if ~exist('conc','var') conc=79; end; %39,79 or empty (1)
if ~exist('probes','var'); probes= 11:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= '19'; end ;  
if ~exist('WaveType','var'); WaveType= 'Regular'; end ;  

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end;
if ~exist('wbin','var'); wbin =3; end ;  

if ~exist('t0description','var'); t0description ='waves reach x'; end ; 
if ~exist('t1description','var')
    if strcmp(WaveType,'Regular')
        if probes(1) == 1
            t1description ='MIZ ref waves reach x';
        else
            t1description ='beach ref waves reach x';
        end
    else
        t1description ='final waves reach x';
    end
end 

if ~exist('OutputStr','var'); OutputStr ='TargetSpectra Mean_AllProbe_Spectra Ind_Spectra_Whole'; end ;  


out_str = ' ''dummy'' '; out_val = ' 0 ';

for i = 1:length(probes)
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probes(i));

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

    %Average of Spectra over time intervals between start time and end time
    Sn1(:,i) = (mean(Sn_mat(:,jj0:jj1), 2));
    
    if strfind(OutputStr,'Ind_Spectra_Whole')
        %Write a matrix for each
        SnA{i} = Sn_mat;
    end
    
end

MeanSpectra  = mean(Sn1,2);
StdSpectra = std(Sn1,[],2);
MinSpectra = min(Sn1,[],2);
MaxSpectra = max(Sn1,[],2);


%Outputs - determine outputs
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
  
  AvgSpecFilt =(1./(period(jj - wbin) - period(jj + wbin))).* trapz(period(jj + [wbin:-1:-wbin]),MeanSpectra(jj + [wbin:-1:-wbin]));
  StdSpecFilt =(1./(period(jj - wbin) - period(jj + wbin))).* trapz(period(jj + [wbin:-1:-wbin]),StdSpectra(jj + [wbin:-1:-wbin]));
  MinSpecFilt = MinSpectra(jj);
  MaxSpecFilt = MaxSpectra(jj);
  
 out_str = [out_str '; ''PerFilt'' '];
 out_val = [out_val '; PerFilt'];
 out_str = [out_str '; ''AvgSpecFilt'' '];
 out_val = [out_val '; AvgSpecFilt'];
 out_str = [out_str '; ''StdSpecFilt'' '];
 out_val = [out_val '; StdSpecFilt'];
 out_str = [out_str '; ''MinSpecFilt'' '];
 out_val = [out_val '; MinSpecFilt'];
 out_str = [out_str '; ''MaxSpecFilt'' '];
 out_val = [out_val '; MaxSpecFilt'];
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

