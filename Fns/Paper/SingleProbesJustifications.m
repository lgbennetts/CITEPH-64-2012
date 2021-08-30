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

function fn_StatisticsOverTimeOneProbe()

if ~exist('conc','var') conc=0; end %39,79 or empty (1)
if ~exist('probe','var');     probe= 1 ;end; 
if ~exist('TestName','var');     TestName= '17' ;end;  %'14' , '6', '16'
if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp,WaveType);'; end

if ~exist('wbin','var');     wbin= 3 ;end; 

close all;

t0description ='waves reach x';
t1description ='final waves reach x';


[ProbeLocXY,tm,disp,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName, probe);
dt = tm(2) - tm(1);

Tp = c_pram.period / 10.0;
%
Hm = c_pram.wave_height /100.0;

Description = ['Concentration ',num2str(conc) '%  , Probe : ', num2str(probe),'  Tp =', num2str(Tp), '[Hz]  Hs =', num2str(Hm), '[m]'  ];

eval(Tpers);
%Tpers = 20;
T_windowSec = round(Tpers*Tp);
T_window = T_windowSec/dt;
SamplingRate = floor(1./dt);
[Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
periods = 1./ff;
FourierWindow_Sec = FourierWindow/SamplingRate;

ffp = ff *Tp;

%Figure for Hs over time;

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType);
t0index = find(strcmp({Tind.description},t0description));
t0 = Tind(t0index).time + T_windowSec;
t1index = find(strcmp({Tind.description},t1description));
t1 = Tind(t1index).time - T_windowSec;

% indextstart=max(find(Trecord<tstart));
jj0=find(tm_vec>t0,1);
jj1=find(tm_vec<t1,1,'last');
t0jj0 = tm_vec(jj0);
t1jj1 = tm_vec(jj1);

MeanSpectra = mean(Sn_mat(:,jj0:jj1), 2);

%S(f) over Time
hmJS =  Hm;
wJS = 2*pi/Tp;
ffJS = 0:0.01:5;
JSspec  = jonswap(2*pi.*ffJS,'wp',wJS,'Hs',hmJS);
ffJSp = ffJS *Tp;
figure();
set(gca,'FontSize',18) 
plot(ffJSp,JSspec,'-k','LineWidth',2, 'DisplayName', 'JonSwap' )
hold on;
plot(ffp,MeanSpectra,'.-b','LineWidth',2, 'DisplayName', 'Mean Spectra', 'MarkerSize',12)
xlim([0,3]);
title(['Average Energy Spectra between t = ',num2str(t0jj0),' [s] ',num2str(t1jj1),' [s] (',num2str(jj1 - jj0),' Time Windows) ', Description]);
xlabel('Frequency / Fp')
ylabel('S(f)')
legend();

%Figure - Hs and Tp
figure()
set(gca,'FontSize',18) 

[~,jjT] = min(abs(Tp - periods));

%Hs
sgtitle(['Key Statistics Over Time ', Description])
subplot(1,2,1);
m0Array = trapz(ff, Sn_mat);
m0Mean = trapz(ff, MeanSpectra);
p1 = plot(tm_vec,4*sqrt(m0Array),'-b', 'DisplayName', 'Hs of Time Windows','LineWidth',2);
hold on;
% p1a = plot(tm_vec,sqrt(2.*sum(Sn_mat(jjT + [-1:1:1],:),1)),'-k', 'DisplayName', 'Peak Height Around','LineWidth',2);
p2 = plot(tm_vec,0*tm_vec + Hm,'--r', 'DisplayName', 'Target Hs','LineWidth',2);
p3 = plot(tm_vec,0*tm_vec + 4*sqrt(m0Mean),'--k', 'DisplayName', 'Mean Spectra Hs','LineWidth',2);
p4 = plot([t0jj0 t0jj0], [0, max(4*sqrt(m0Array))], ':k' , 'DisplayName', ['Time Interval Boundaries'],'LineWidth',2);
plot([t1jj1 t1jj1], [0, max(4*sqrt(m0Array))], ':k','LineWidth',2)
xlabel('Mean Time [s]')
ylabel('Hs [m]')
legend([p1 p2 p3 p4])
%Tp
subplot(1,2,2);
[~,TpArrayI] = max(Sn_mat,[],1);
Tps = periods(TpArrayI);
[~,TpMSI] = max(MeanSpectra);
TpMS = periods(TpMSI);
p1 = plot(tm_vec,Tps,'-b', 'DisplayName', 'Tp of Time Windows','LineWidth',2);
hold on;
p2 = plot(tm_vec,0*tm_vec + Tp,'--r', 'DisplayName', 'Target Tp','LineWidth',2);
p3 = plot(tm_vec,0*tm_vec + TpMS,'--k', 'DisplayName', 'Mean Spectra Tp','LineWidth',2);
p4 = plot([t0jj0 t0jj0], [min(Tps), max(Tps)], ':k' , 'DisplayName', 'Time Interval Boundaries','LineWidth',2);
plot([t1jj1 t1jj1], [min(Tps), max(Tps)], ':k','LineWidth',2)
xlabel('Mean Time [s]')
ylabel('Tp [s]')
legend([p1 p2 p3 p4])

Hserror = abs(4*sqrt(m0Mean) - Hm) / Hm
Tperror = abs(TpMS - Tp)/ Tp


%Raw Signal
% figure();
% plot(tm,disp,'-k')
% hold on;
% p1 = plot([t0jj0 t0jj0], [min(disp),max(disp)], '--b' , 'DisplayName', ['Time Interval Begin']);
% p2 = plot([t1jj1 t1jj1], [min(disp),max(disp)], '--r' , 'DisplayName', ['Time Interval End']);
% hold off;
% title(['Raw Signal ' , Description]);
% xlabel('Time [s]')
% ylabel('Displacement [m]')
% legend([p1,p2])

% % Significant Wave Height over Time
% figure();
% m0Array = trapz(ff, Sn_mat);
% m0Mean = trapz(ff, MeanSpectra);
% plot(tm_vec,4*sqrt(m0Array), 'DisplayName', 'From Spectra');
% hold on;
% plot(tm_vec,0*tm_vec + Hm,'--k', 'DisplayName', 'Experimental Target');
% plot(tm_vec,0*tm_vec + 4*sqrt(m0Mean),'--r', 'DisplayName', 'Mean Spectra Hs');
% plot([t0jj0 t0jj0], [0, max(4*sqrt(m0Array))], ':b' , 'DisplayName', ['Time Interval Begin'])
% plot([t1jj1 t1jj1], [0, max(4*sqrt(m0Array))], ':r' , 'DisplayName', ['Time Interval End'])
% xlabel('average time (s)')
% ylabel('Hs = $4 \sqrt{m_0}$ [m]','Interpreter','latex')
% title(['Hs Over Time ', Description])
% legend()
% 
% 
% % Tp over Time
% figure();
% m1Array = trapz(ff, ff.*Sn_mat);
% m1Mean = trapz(ff, ff.*MeanSpectra);
% 
% m01Array = m0Array ./ m1Array;
% m01Mean =  m0Mean./m1Mean;
% 
% JSm0 = (Hm/4)^2;
% JSm1 = trapz(ffJS, ffJS.*JSspec);
% JSTm01 = JSm0 ./ JSm1;
% plot(tm_vec,m01Array, 'DisplayName', 'From Spectra');
% hold on;
% plot(tm_vec,0*tm_vec + JSTm01,'--k', 'DisplayName', 'Experimental Target');
% plot(tm_vec,0*tm_vec + m01Mean,'--r', 'DisplayName', 'Mean Tp of All Spectra In Steady State');
% plot([t0jj0 t0jj0], [0, 1.5], ':b' , 'DisplayName', ['Time Interval Begin'])
% plot([t1jj1 t1jj1], [0,1.5], ':r' , 'DisplayName', ['Time Interval End'])
% xlabel('average time (s)')
% ylabel('Tm01')
% title(['Tm01 Over Time ', Description])
% legend()

%Figure - m0,m1,m2,m3
% m2Array = trapz(ff, ff.^2.*Sn_mat);
% m2Mean = trapz(ff, ff.^2.*MeanSpectra);
% JSm2 = trapz(ffJS, ffJS.^2.*JSspec);


% figure();
% sgtitle(['Moments Over Time ', Description])
% subplot(3,1,1);
% p1 = plot(tm_vec,m0Array,'-b', 'DisplayName', 'm0');
% hold on;
% p2 = plot(tm_vec,0*tm_vec + m0Mean,'--b', 'DisplayName', 'Mean m0 of All Spectra In Steady State');
% p3 = plot(tm_vec,0*tm_vec + JSm0,'--r', 'DisplayName', 'm0 of JonSwap');
% p4 = plot([t0jj0 t0jj0], [0, max(m0Array)], '--k' , 'DisplayName', 'Time Interval Begin');
% plot([t1jj1 t1jj1], [0, max(m0Array)], '--k', 'DisplayName','')
% xlabel('Average Time (s)')
% ylabel('Zero Moment')
% legend([p1 p2 p3 p4])
% subplot(3,1,2);
% p1 = plot(tm_vec,m1Array,'-b', 'DisplayName', 'm1');
% hold on;
% p2 = plot(tm_vec,0*tm_vec + m1Mean,'--b', 'DisplayName', 'Mean m1 of All Spectra In Steady State');
% p3 = plot(tm_vec,0*tm_vec + JSm1,'--r', 'DisplayName', 'm1 of JonSwap');
% p4 = plot([t0jj0 t0jj0], [0, max(m1Array)], '--k' , 'DisplayName', 'Time Interval Begin');
% p5 = plot([t1jj1 t1jj1], [0, max(m1Array)], '--k', 'DisplayName','');
% xlabel('Average Time (s)')
% ylabel('First Moment')
% legend([p1 p2 p3 p4])
% subplot(3,1,3);
% p1 = plot(tm_vec,m2Array,'-b', 'DisplayName', 'm2');
% hold on;
% p2 = plot(tm_vec,0*tm_vec + m2Mean,'--b', 'DisplayName', 'Mean m2 of All Spectra In Steady State');
% p3 = plot(tm_vec,0*tm_vec + JSm2,'--r', 'DisplayName', 'm2 of JonSwap');
% p4 = plot([t0jj0 t0jj0], [0, max(m2Array)], '--k' , 'DisplayName', 'Time Interval Begin');
% p5 = plot([t1jj1 t1jj1], [0, max(m2Array)], '--k', 'DisplayName','');
% xlabel('Average Time (s)')
% ylabel('Second Moment')
% legend([p1 p2 p3 p4])

return


function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

