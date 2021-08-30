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

function fn_CompareProbesToInitial(conc, probe,TestName)


%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Which user?
[user_number,user_name] = citeph_get_user();

if user_number == 0
    OTHER_USR=1;
end


%% GENERAL
%pick expriment
if ~exist('conc','var') conc=39; end %39,79 or empty (1)
if ~exist('probe','var');     probe= 17 ;end; 
if ~exist('TestName','var');     TestName= '1' ;end; 

if ~exist('wbin','var');   wbin =3; end


if ~exist('data_out','var')
 data_out.name={'amp-harmo-steady-1st'};
 data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
end

if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

%Read Data in Calib 1, at probe
[ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(1,TestName, probe);

if Success == 0
   [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, 1);
end

close all;
Tp = c_pram.period/10.0;
dt = tm(2) - tm(1);

description =  [WaveType ' waves;' ...
    ' Concentration =' num2str(conc) ';' ...
   ' Hs=' num2str(10*c_pram.wave_height) '[mm];' ...
   ' Tm=' num2str(c_pram.period/10) ' [s];' ...
   ' probe(s)=' int2str(probe)] ;


[Energy_Baseline,Peaks_Baseline,Regions_Baseline,TgTps_Baseline] = fn_ExtractAllPeaksProbes(tm,disp,Tp,conc,ProbeLocXY(1),dt);
clear ProbeLocXY tm disp c_pram WaveType;



% % %Read Data in Probe
% [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probe);
% % 
% [Energy_Pbx] = fn_EnergyAtPeaks(tm,disp,conc,ProbeLocXY(1),dt,Peaks_Baseline,Regions_Baseline,TgTps_Baseline);


% % Tind = fn_Tind(conc,Tp,ProbeLocXY(1)); %Slower
% Tind = fn_Tind(conc,Tp,ProbeLocXY(1))
% % TindB = fn_Tind(conc,max(Peaks_Pb1),ProbeLocXY(1)); % Faster
% 
% %pick t0 as slowest possible time
% %pick t1 as fastest posssible time
% 
% t0 = Tind(7).time;
% t1 = Tind(8).time;
% 
% 
% dtlength = (t1 - t0) / dt;
% 
% timewindowlength = 2^floor(log2(dtlength));
% SI = floor((t0 -tm(1)) / dt);
% EI = SI + timewindowlength-1;
% 
% %% Fourier Transform
% [Sn,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
% periods = 1./ff;
% 
% figure();
% plot(periods,Sn,'-');
% hold on;
% for i = 1: length(Peaks_Baseline)
%     plot([Peaks_Baseline(i),Peaks_Baseline(i)],[0,max(Sn)],'--k');
%     txt = ['T = ' num2str(Peaks_Baseline(i))];
%     text(Peaks_Baseline(i)*1.01,0.5*max(Sn),txt)
% end
% title([description '  Found Spectral Peaks (From Probe 1) and Spectrogram for Ideal Time Window For Tp']);
% xlabel('Period (s)');
% xlabel('Spectogram ');
% hold off;

%Are these good peaks?
%plot against FFT
% 
% Energy_Pbx ./ Energy_Baseline

%Other Probe

return


function Tind = fn_Tind(conc,Tp,X)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return

