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

function fn_InvestigateSpectra(conc, probe,TestName)


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
if ~exist('conc','var') conc=1; end %39,79 or empty (1)
if ~exist('probe','var');     probe= 1 ;end; 
if ~exist('TestName','var');     TestName= '15' ;end; 

if ~exist('wbin','var');   wbin =3; end


if ~exist('data_out','var')
 data_out.name={'amp-harmo-steady-1st'};
 data_out.tint='[t0,t1] = fn_tint(tvec,t0,Tp,Tpers);';
end

if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

%Read Data
[ProbeLocXY,tm,disp,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName, probe);

Tp = c_pram.period / 10.0;

description =  [WaveType ' waves;' ...
   ' Hs=' num2str(10*c_pram.wave_height) '[mm];' ...
   ' Tm=' num2str(c_pram.period/10) ' [s];' ...
   ' probe(s)=' int2str(probe)] ;

%Plot
% figure();
% plot(tm,disp);
% title([description 'Signal' ])
% xlabel('Time(s)')
% ylabel('Displacement (mm)')

%% Fourier Transform

%
tmRed=find(and(tm<=250, tm >= 50));

[Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,Tpers);


periods = 1./ff;

 %%%: Determine time window:
 

%find mean of the peaks throughout time interval
% [maxvalue, PeakPeriod] = max(Sn_mat,[],1);
% TpC = mean(periods(PeakPeriod));

%Find closest freuqency to the target - want to plot it, +-1, +-2, +-3
[~,jj]=min(abs(periods-Tp));
TpC = periods(jj);

TotalWaveEnergyAtTargPer = sum(Sn_mat(jj+[-1:1],:));

[t0,t1] = fn_FindStationaryWindow(tm_vec,TotalWaveEnergyAtTargPer,Tp,ProbeLocXY(1),conc)


%TpC = Tp;
 
% Tind = fn_Tind(conc,TpC,ProbeLocXY(1),tm);
% tvec=[]; count=1;
% for loop=1:length(Tind)
%     if strfind(Tind(loop).description,'x')
%         tvec(count)=Tind(loop).time;
%         tvecdes{count}=Tind(loop).description; count=count+1;
%     end
% end % end loop Tind
% clear count Tind

% [t0,it0]=min(tvec); tvec(it0)=[]; tvecdes(it0)=[]; clear it0
% if strcmp(data_out.tint,'default')
%     t1=min(tvec); clear tvec
% else
%     eval(data_out.tint)
%     clear tvec
% end % end if tint='default'
% tt=find(and(tm_vec>=t0,tm_vec<=t1));



% maxSnMat = max(max(Sn_mat(:,:)));
% 
% figure();
% plot(ff(:),Sn_mat(:,tt));
% title([description 'Wave Spectrum At Mean Times All Time'])
% xlabel('Frequency (Hz)')
% ylabel('Wave Spectrum E(t_m)')
% hold on;
% plot([ff(jj) , ff(jj) ],[0,1.1*maxSnMat],'--k');
% 
% plot([ff(jj-1) , ff(jj-1) ],[0,1.1*maxSnMat],'--b');
% plot([ff(jj+1) , ff(jj+1) ],[0,1.1*maxSnMat],'--b');
% 
% plot([ff(jj-3) , ff(jj-3) ],[0,1.1*maxSnMat],'--r');
% plot([ff(jj+3) , ff(jj+3) ],[0,1.1*maxSnMat],'--r');
% hold off;
% 
% [MaxValue,PeakPeriods] = maxk(Sn_mat,5,1);
% 
% size(tm_vec)

% figure();
% plot(tm_vec(tt),ff(PeakPeriods(1,tt)),'-b');
% hold on;
% % plot(tm_vec,ff(PeakPeriods(2,:)),'-r');
% % plot(tm_vec,ff(PeakPeriods(3,:)),'-k');
% % plot(tm_vec,ff(PeakPeriods(4,:)),'--g');
% % plot(tm_vec,ff(PeakPeriods(5,:)),'--k');
% title([description 'Peak Wave Spectra']);
% xlabel('Mean Time (s)')
% ylabel('Wave Spectrum E(t_m)')

%Amplitude +-1 away from target peak
% figure();
% TotalWaveEnergyAtTargPer = sum(Sn_mat(jj+[-1:1],:));
% plot(tm_vec(:,1),TotalWaveEnergyAtTargPer(:),'-b');
% hold on;
% plot([t0 , t0 ],[0,1.1*max(TotalWaveEnergyAtTargPer)],'--k');
% plot([t1 , t1 ],[0,1.1*max(TotalWaveEnergyAtTargPer)],'--k');
%More focused window
 
% figure();
% plot(ff(:),Sn_mat(:,tt));
% title([description 'Wave Spectrum At Between' int2str(t0) '  and  ' int2str(t1)])
% xlabel('Freuqency (Hz)')
% ylabel('Wave Spectrum E(t_m)')


% figure();
% [maxvalue, PeakPeriod] = max(Sn_mat,[],1);
% plot(tm_vec,periods(PeakPeriod));
% title([description 'Peak Wave Spectra']);
% xlabel('Mean Time (s)')
% ylabel('Wave Spectrum E(t_m)')
return


function Tind = fn_Tind(conc,Tp,X,tm)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn');
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration');
end

 return

