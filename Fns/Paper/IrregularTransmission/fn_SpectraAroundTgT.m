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

function out = fn_SpectraAroundTgT(Tp,TpLB,TpUB,conc,TestName,probes,OutputStr,TpersStr,t0description,t1description)

% 
if ~exist('conc','var') conc=79; end; %39,79 or empty (1)
if ~exist('probes','var'); probes= 11:20 ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= '19'; end ;  
if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end;
if ~exist('t0description','var'); t0description ='waves reach x'; end ;  
if ~exist('t1description','var'); t1description ='final waves reach x'; end ;  

if ~exist('OutputStr','var'); OutputStr ='TargetSpectra Mean_AllProbe_Spectra Ind_Spectra_Whole'; end ;  

close all;

out_str = ' ''dummy'' '; out_val = ' 0 ';

for i = 1:length(probes)
 %Read Data in Calib 1, at probe
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probes(i));
    
    % Tp = c_pram.period/10.0;
    if ~exist('Tp','var'); Tp = c_pram.period/10.0; end ;  
    if ~exist('TpLB','var'); TpLB = Tp/2; end ;
    if ~exist('TpUB','var'); TpUB = 2*Tp; end ;
    
    dt = tm(2) - tm(1); 

    eval(TpersStr);
    T_windowSec = round(Tpers*Tp);
    SamplingRate = floor(1./dt);
    T_window = SamplingRate*T_windowSec;

    Hs = c_pram.wave_height/100.0;
    
%     Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType); %Slower
    TindLB = fn_Tind(conc,TpLB,ProbeLocXY(1),WaveType); %Slower
    TindUB = fn_Tind(conc,TpUB,ProbeLocXY(1),WaveType); %Slower
    
    t0index = find(strcmp({TindLB.description},t0description));
    t0 = TindLB(t0index).time + T_windowSec;
    t1index = find(strcmp({TindUB.description},t1description));
    t1 = TindUB(t1index).time - T_windowSec ;

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
 out_str = [out_str '; ''MeanSpectra'' '];
 out_val = [out_val '; MeanSpectra'];
 out_str = [out_str '; ''Period'' '];
 out_val = [out_val '; period'];
end


if strfind(OutputStr,'TargetSpectra')
 PerFilt = period(and(period > TpLB, period < TpUB));
 SpecFilt = MeanSpectra(and(period > TpLB, period < TpUB));
 
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

