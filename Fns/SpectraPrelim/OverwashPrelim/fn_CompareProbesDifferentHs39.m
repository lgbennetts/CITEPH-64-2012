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

function fn_StatisticsOverTime()

if ~exist('conc','var') conc =39; end %39,79 or empty (1)
if ~exist('probeT','var');     probeT= 11 ;end; 
if ~exist('probeR','var');     probeR= 9 ;end; 
if ~exist('TestName','var');     TestName= {'2','7'} ;end;  %'2' and '7' %TestName= {'3','6'}
if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

if ~exist('wbin','var');     wbin= 3 ;end; 

close all;

%First Test

%Reflection Probe

[ProbeLocXY,tm1R,disp1R,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName{1}, probeR);

Tp(1) = c_pram.period / 10.0;
Hm(1) = c_pram.wave_height /100.0;
eval(Tpers);
Tpers = 100;
T_window = Tpers*Tp(1);

[Sn_mat1R,tm_vec1R,ff1R] = fn_JustFourierMovingWindow(tm1R,disp1R,Tp(1),T_window);
period1R = 1./ff1R;

%Transmission Probe

[ProbeLocXY,tm1T,disp1T,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName{1}, probeT);
[Sn_mat1T,tm_vec1T,ff1T] = fn_JustFourierMovingWindow(tm1T,disp1T,Tp(1),T_window);
period1T = 1./ff1T;


Test1Description = ['Conc ',num2str(conc) '%  , ProbeR : ', num2str(probeR),' ProbeT : ', num2str(probeT),'  Tp =', num2str(Tp(1)), '[Hz]  Hs =', num2str(Hm(1)), '[m]'  ];


%Second Test
%Reflection Probe
[ProbeLocXY,tm2R,disp2R,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName{2}, probeR);

Tp(2) = c_pram.period / 10.0;
Hm(2) = c_pram.wave_height /100.0;
T_window = Tpers*Tp(2);

[Sn_mat2R,tm_vec2R,ff2R] = fn_JustFourierMovingWindow(tm2R,disp2R,Tp(2),T_window);
period2R = 1./ff2R;

%Transmission Probe

[ProbeLocXY,tm2T,disp2T,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName{2}, probeT);
[Sn_mat2T,tm_vec2T,ff2T] = fn_JustFourierMovingWindow(tm2T,disp2T,Tp(2),T_window);
period2T = 1./ff2T;

Test2Description = ['Conc ',num2str(conc) '%  , ProbeR : ', num2str(probeR),' ProbeT : ', num2str(probeT),'  Tp =', num2str(Tp(2)), '[Hz]  Hs =', num2str(Hm(2)), '[m]'  ];

Description = ['Conc ',num2str(conc) '%  , ProbeR : ', num2str(probeR),' ProbeT : ', num2str(probeT),'  Tp =', num2str(Tp), '[Hz]  Hs =', num2str(Hm), '[m]'  ];


%Raw Signal - Test 1
figure();
plot(tm1R,disp1R,'-r','DisplayName','Test 1 Reflection Probe')
hold on;
plot(tm1T,disp1T,'-k','DisplayName','Test 1 Transmission Probe')
hold off;
title(['Raw Signal ' , Test1Description]);
xlabel('time (s)')
ylabel('displacement (m)')
legend();

%Raw Signal
figure();
plot(tm2R,disp2R,'-r','DisplayName','Test 2 Reflection Probe')
hold on;
plot(tm2T,disp2T,'-k','DisplayName','Test 2 Transmission Probe')
hold off;
title(['Raw Signal ' , Test2Description]);
xlabel('time (s)')
ylabel('displacement (m)')
legend();

%Spectra - 
targtime = 100;

[~,jj]=min(abs(tm_vec1R-targtime));
Sn_1R = Sn_mat1R(:,jj);

[~,jj]=min(abs(tm_vec1T-targtime));
Sn_1T = Sn_mat1T(:,jj);

[~,jj]=min(abs(tm_vec2R-targtime));
Sn_2R = Sn_mat2R(:,jj);

[~,jj]=min(abs(tm_vec2T-targtime));
Sn_2T = Sn_mat2T(:,jj);

figure();
plot(ff1R*Tp(1),Sn_1R,'-r','DisplayName','Test 1 Reflection Probe')
hold on;
plot(ff1T*Tp(1),Sn_1T,'-k','DisplayName','Test 1 Transmission Probe')
plot(ff2R*Tp(1),Sn_2R,'-b','DisplayName','Test 2 Reflection Probe')
plot(ff2T*Tp(1),Sn_2T,'-g','DisplayName','Test 2 Transmission Probe')
hold off;
set(gca, 'YScale', 'log')
title(['Spectra ' , Description, ' Time Window ', num2str(T_window) ,' [s]']);
xlabel('frequency / inc frequency')
ylabel('Spectra')
% xlim([0,2*Tp(1)])
legend();

figure();
plot(period1R,Sn_1R,'-r','DisplayName','Test 1 Reflection Probe')
hold on;
plot(period1T,Sn_1T,'-k','DisplayName','Test 1 Transmission Probe')
plot(period2R,Sn_2R,'-b','DisplayName','Test 2 Reflection Probe')
plot(period2T,Sn_2T,'-g','DisplayName','Test 2 Transmission Probe')
hold off;
title(['Spectra ' , Description, ' Time Window ', num2str(T_window) ,' [s]']);
xlabel('period')
ylabel('Spectra')
 xlim([0,2*Tp(1)])
legend();

% figure();
% 
% Sn_1TMax = max(Sn_1T);
% Sn_2TMax = max(Sn_2T);
% PerC = 0.01;
% 
% plot(period1R(Sn_1T > PerC*Sn_1TMax),Sn_1T(Sn_1T > PerC*Sn_1TMax)./Sn_1R(Sn_1T > PerC*Sn_1TMax),'.r','DisplayName','Test 1 Transmission Ratio')
% hold on;
% plot(period2R(Sn_2T > PerC*Sn_2TMax),Sn_2T(Sn_2T > PerC*Sn_2TMax)./Sn_2R(Sn_2T > PerC*Sn_2TMax),'.b','DisplayName','Test 2 Transmission Ratio')
% hold off;
% title(['Transmission Ratio ' , Description, ' Time Window ', num2str(T_window) ,' [s]']);
% xlabel('period (s)')
% ylabel('Transmission Ratio')
% xlim([0,2*Tp(1)])
% legend();


%Spectra Over Time and Power in
% figure();
% wbin = 5;
% [~,jj1R]=min(abs(period1R-Tp));
% jj1R = jj1R(1);
% 
% [~,jj1T]=min(abs(period1T-Tp));
% jj1T = jj1T(1);
% 
% [~,jj2R]=min(abs(period2R-Tp));
% jj2R = jj2R(1);
% 
% [~,jj2T]=min(abs(period2T-Tp));
% jj2T = jj2T(1);
% 
% plot(tm_vec1R(:),sum(Sn_mat1R(jj1R + [-wbin:wbin],:),1),'-b','DisplayName','Test 1 Reflection Probe')
% hold on;
% plot(tm_vec1T(:),sum(Sn_mat1T(jj1T + [-wbin:wbin],:),1),'-r','DisplayName','Test 1 Transmission Probe')
% plot(tm_vec2R(:),sum(Sn_mat2R(jj2R + [-wbin:wbin],:),1),'-k','DisplayName','Test 2 Reflection Probe')
% plot(tm_vec2T(:),sum(Sn_mat2T(jj2T + [-wbin:wbin],:),1),'-g','DisplayName','Test 2 Transmission Probe')
% 
% title('Amplitude At Peak Frequency (Sum Neighbours)');
% xlabel('average time (s)')
% ylabel('Energy')
% legend();

return



