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
LBFac = 0.75;
UBFac = 1.5;


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


%Figure for Hs over time;

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType);
TindLB = fn_Tind(conc,LBFac*Tp,ProbeLocXY(1),WaveType);
TindUB = fn_Tind(conc,UBFac*Tp,ProbeLocXY(1),WaveType);
t0index = find(strcmp({TindLB.description},t0description));
t0 = TindLB(t0index).time + FourierWindow_Sec;
t1index = find(strcmp({TindUB.description},t1description));
t1 = TindUB(t1index).time - FourierWindow_Sec;


[~,jj0]=min(abs(tm_vec-t0));
[~,jj1]=min(abs(tm_vec-t1));


t0jj0 = tm_vec(jj0);
t1jj1 = tm_vec(jj1);

MeanSpectra = mean(Sn_mat(:,jj0:jj1), 2);

% m = 1;
% 
% MeanSpectra = downsample(movmean(MeanSpectra, m),m);
% ffds = downsample(ff,m);

%Raw Signal
figure();
plot(tm,disp,'-k')
hold on;
p1 = plot([t0jj0 t0jj0], [min(disp),max(disp)], '--b' , 'DisplayName', ['Time Interval Begin']);
p2 = plot([t1jj1 t1jj1], [min(disp),max(disp)], '--r' , 'DisplayName', ['Time Interval End']);
hold off;
title(['Raw Signal ' , Description]);
xlabel('Time [s]')
ylabel('Displacement [m]')
legend([p1,p2])

%Surf Map of Signal and Amplitude Over Time
%Regular
figure();
sgtitle(['Amplitude ' ,Description, '  Time Window = ', num2str(FourierWindow_Sec),' [s]']);
h1=subplot(2,1,1);
surf(h1,tm_vec(:),ff*Tp,log10(sqrt(2*Sn_mat)))
shading interp
view(h1,[0,90])
xlim(h1,[tm_vec(1), tm_vec(end)])
ylim(h1,[0 5]);
cb=colorbar; set(get(cb,'ylabel'),'String','log_{10}(A)'); 
cb=get(cb,'ylabel'); set(cb,'fontsize',16); clear cb
xlabel('mean time [s]')
ylabel('frequency / fp')
title('Spectra for All Time Windows')

h2=subplot(2,1,2);
plot(h2,ff*Tp,sqrt(2*MeanSpectra),'-b','LineWidth',2, 'DisplayName', 'Data')
hold on;
wJS = 2*pi/Tp;
hmJS =  Hm;
ffJS = 0:0.01:5;
JSspec  = jonswap(2*pi.*ffJS,'wp',wJS,'Hs',hmJS);
JSspec = sqrt(2*JSspec);

plot(h2,ffJS*Tp,JSspec,'--k','LineWidth',2, 'DisplayName', ['JonSwap wp: ', num2str(wJS), '[Hz]' , ' Hs ($4\sqrt{m_0}$) ', num2str(round(hmJS,4)), '[m]' ])

hold off;
xlim(h2,[0 5]);
title(['Average Amplitude Spectra between t = ',num2str(tm_vec(jj0)),' [s] ',num2str(tm_vec(jj1)),' [s] (',num2str(jj1 - jj0),' Time Windows)']);
xlabel('frequency / fp')
ylabel('amplitude (m)')
legend('Interpreter','latex');


%Surf Map of Signal and Energy Over Time
%Regular
figure();
sgtitle(['Energy ' ,Description, '  Time Window = ', num2str(FourierWindow_Sec),' [s]']);
h1=subplot(2,1,1);
surf(h1,tm_vec(:),ff*Tp,log10(Sn_mat))
shading interp
view(h1,[0,90])
xlim(h1,[tm_vec(1), tm_vec(end)])
ylim(h1,[0 5]);
cb=colorbar; set(get(cb,'ylabel'),'String','log_{10}(S)'); 
cb=get(cb,'ylabel'); set(cb,'fontsize',16); clear cb
xlabel('average time (s)')
ylabel('frequency / fp')
title('Spectra for All Time Windows')

h2=subplot(2,1,2);
plot(h2,ff*Tp,MeanSpectra,'-b','LineWidth',2, 'DisplayName', 'Data')
hold on;
hmJS =  Hm;
ffJS = 0:0.01:5;
JSspec  = jonswap(2*pi.*ffJS,'wp',wJS,'Hs',hmJS);
plot(h2,ffJS*Tp,JSspec,'--k','LineWidth',2, 'DisplayName', ['JonSwap wp: ', num2str(wJS), '[Hz]' , '  Hs ($4\sqrt{m_0}$) ', num2str(round(hmJS,4)), '[m]' ])
hold off;
xlim(h2,[0 ,5]);
title(['Average Energy Spectra between t = ',num2str(t0jj0),' [s] ',num2str(t1jj1),' [s] (',num2str(jj1 - jj0),' Time Windows)']);
xlabel('frequency / fp')
ylabel('Spectra')
legend('Interpreter','latex');

% Significant Wave Height over Time
figure();
m0Array = trapz(ff, Sn_mat);
m0Mean = trapz(ff, MeanSpectra);
plot(tm_vec,4*sqrt(m0Array), 'DisplayName', 'From Spectra');
hold on;
plot(tm_vec,0*tm_vec + Hm,'--k', 'DisplayName', 'Experimental Target');
plot(tm_vec,0*tm_vec + 4*sqrt(m0Mean),'--r', 'DisplayName', 'Mean Spectra Hs');
plot([t0jj0 t0jj0], [0, max(4*sqrt(m0Array))], ':b' , 'DisplayName', ['Time Interval Begin'])
plot([t1jj1 t1jj1], [0, max(4*sqrt(m0Array))], ':r' , 'DisplayName', ['Time Interval End'])
xlabel('average time (s)')
ylabel('Hs = $4 \sqrt{m_0}$ [m]','Interpreter','latex')
title(['Hs Over Time ', Description])
legend()



return


function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

