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

if ~exist('conc','var') conc=39; end %39,79 or empty (1)
if ~exist('probe','var');     probe= 1 ;end; 
if ~exist('TestName','var');     TestName= '9' ;end;  %'14' , '6'
if ~exist('Tpers','var');  Tpers='Tpers=fn_Tpers(Tp);'; end

if ~exist('wbin','var');     wbin= 3 ;end; 

close all;

t0description ='waves reach x';
t1description ='final waves reach x';



[ProbeLocXY,tm,disp,c_pram,WaveType] = fn_FindAndReadProbe(conc,TestName, probe);

Tp = c_pram.period / 10.0;
%
Hm = c_pram.wave_height /100.0;

Description = ['Concentration ',num2str(conc) '%  , Probe : ', num2str(probe),'  Tp =', num2str(Tp), '[Hz]  Hs =', num2str(Hm), '[m]'  ];

%eval(Tpers);
Tpers = 50;
T_window = Tpers*Tp;

[Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
periods = 1./ff;

%Raw Signal
figure();
plot(tm,disp,'-k')
title(['Raw Signal ' , Description]);
xlabel('time (s)')
ylabel('displacement (m)')

%Surf Map of Signal and Amplitude Over Time
%Regular
figure();
sgtitle([Description, '  Time Window = ', num2str(T_window),' [s]']);
h1=subplot(2,1,1);
surf(h1,tm_vec(:),periods,log10(sqrt(2*Sn_mat)))
shading interp
view(h1,[0,90])
xlim(h1,[tm_vec(1), tm_vec(end)])
ylim(h1,[0 2*Tp]);
cb=colorbar; set(get(cb,'ylabel'),'String','log_{10}(A)'); 
cb=get(cb,'ylabel'); set(cb,'fontsize',16); clear cb
xlabel('average time (s)')
ylabel('period (s)')
title('Spectra for All Time Windows')

%Time for
LBFac = 1;
UBFac = 1;

Tind = fn_Tind(conc,Tp,ProbeLocXY(1),WaveType);
TindLB = fn_Tind(conc,LBFac*Tp,ProbeLocXY(1),WaveType);
TindUB = fn_Tind(conc,UBFac*Tp,ProbeLocXY(1),WaveType);
t0index = find(strcmp({TindLB.description},t0description));
t0 = TindLB(t0index).time + T_window;
t1index = find(strcmp({TindUB.description},t1description));
t1 = TindUB(t1index).time - T_window;

[~,jj0]=min(abs(tm_vec-t0));
[~,jj1]=min(abs(tm_vec-t1));

h2=subplot(2,1,2);
MeanSpectra = mean(Sn_mat(:,jj0:jj1), 2);
plot(h2,periods,sqrt(2*MeanSpectra),'-b', 'DisplayName', 'Data')
hold on;

%Area under curve
Cff = ff(and(periods > Tp/2, periods < 2*Tp));
CMeanSpectra = MeanSpectra(and(periods > Tp/2, periods < 2*Tp));
m0 = trapz(Cff, CMeanSpectra);

%Why this factor of 10 difference???????
%What is Hs in jonswap?
%What is Hs in experiment?
wJS = 2*pi/Tp;
hmJS = 4*sqrt(m0); %Hm;
% hmJS = Hm;
JSspec  = jonswap(2*pi./periods,'wp',wJS,'Hs',hmJS);
JSspec = sqrt(2*JSspec);

plot(h2,periods,JSspec,'--k', 'DisplayName', ['JonSwap wp: ', num2str(wJS), '[Hz]' , '  Hs (4sqrt(m0)) ', num2str(round(hmJS,4)), '[m]' ])

hold off;
xlim(h2,[0 2*Tp]);
title(['Average Amplitude Spectra between t = ',num2str(tm_vec(jj0)),' [s] ',num2str(tm_vec(jj1)),' [s]']);
xlabel('average time (s)')
ylabel('amplitude (m)')
legend();

% figure();
% [~,jj]=min(abs(periods-Tp));
% plot(tm_vec(:),sqrt(2*sum(Sn_mat(jj+[-wbin:wbin],:),1)),'-')
% title('Amplitude At Peak Frequency (Sum Neighbours)');
% xlabel('average time (s)')
% ylabel('elevation (m)')

% 
% %Bin Spectra, and a look at energy over time
% StepSize = 0.05;
% 
% BinMids = 0.6:StepSize:1.2;
% 
% BinEdges(1) = -1;
% BinEdges(2:length(BinMids)+1) = BinMids - StepSize/2;
% BinEdges(end + 1) = BinMids(end) + StepSize/2;
% BinEdges(end+1) = 20000;
% 
% [N,edges,bin] = histcounts(periods,BinEdges);
% 
% % BinLarge = repmat(bin,1, size(Sn_mat,2));
% [xx, yy] = ndgrid(bin,1:size(Sn_mat,2));
% SnBin = accumarray([xx(:) yy(:)],Sn_mat(:));
% % C=accumarray([xx(:) yy(:)],B(:));
% figure();
% hold on;
% for i = 2:size(SnBin,1)-1
%     Label = ['Range : ', num2str(BinEdges(i)), ' , ' num2str(BinEdges(i+1)) ];
%     plot(tm_vec(:),(SnBin(i,:)),'DisplayName',Label,'LineWidth',2);
% end
% hold off;
% title(['Energy In Period Ranges ' ,Description , '  Time Window = ', num2str(T_window),' [s]'])
% legend()
% xlabel('average time (s)')
% ylabel('Energy')

return


function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

