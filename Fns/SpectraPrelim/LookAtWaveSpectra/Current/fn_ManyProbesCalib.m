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

function fn_ManyProbesCalib(conc,TestName,Tp,probes)

%Compare calibration to disk
if ~exist('conc','var') conc={0,0}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'17','17'}}; end ;  % 9, 14 %14 %2

% if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  

% if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  % 9, 14


if ~exist('Cols','var'); Cols= {'xr', '^b','sk'}; end ;

if ~exist('wbinwidth','var');  wbinwidth=6; end;
if ~exist('UBfactor','var');  UBfactor=1.5; end;
if ~exist('LBfactor','var');  LBfactor=0.75; end;

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end
if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;

j = 1;
for j = 1: length(TestName)
    

    t0description ='waves reach x';
    t1description ='final waves reach x';

    for i = 1:length(probes{1})
     %Read Data in Calib 1, at probe
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, probes{1}(i));

        dt = tm(2) - tm(1);
        
        Tp = c_pram.period/10.0;
        
        eval(TpersStr);
%         Tpers = 20;
        
        %Tpers = fn_Tpers(Tp);
        T_windowSec = round(Tpers*Tp);
        SamplingRate = floor(1./dt);
        T_window = SamplingRate*T_windowSec;
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs = c_pram.wave_height/100.0;

        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        TindLB = fn_Tind(conc{1},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{1},TpUB,ProbeLocXY(1),WaveType);


        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time + T_windowSec;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time - T_windowSec ;

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate); % fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        t0jj0 = tm_vec(jj0);
        t1jj1 = tm_vec(jj1);
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        
        m0Array = trapz(ff, Sn_mat);
        Hs1(:,i) = 4*sqrt(m0Array);
        Tm1(:,i) = tm_vec - t0;    
        TmE1(1,i) = 0;
        TmE1(2,i) = t1 - t0;

        
        %Raw signal figure
        if i == 1
            figure();
            plot(tm,disp,'-k')
            hold on;
            p1 = plot([t0jj0 t0jj0], [min(disp), max(disp)], '--b' , 'DisplayName', ['Time Interval Begin']);
            p2 = plot([t1jj1 t1jj1], [min(disp), max(disp)], '--r' , 'DisplayName', ['Time Interval End']);
            legend([p1 p2])
            ylabel('Displacement [m]')
            xlabel('Time [s]')
            title('Raw Signal')
            hold off;
        end

        %% Fourier Transform

        ff1{i} = ff;

        SnS1{i} = Sn1;
        Times1(i,:) = [tm_vec(jj0)- t0;tm_vec(jj1)- t0];
        Times1US(i,:) = [tm_vec(jj0);tm_vec(jj1)];

        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{j}{2}, probes{2}(i));

        Tp = c_pram.period/10.0;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{2},Tp,ProbeLocXY(1),WaveType); %Slower
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs = c_pram.wave_height/100.0;
        
        TindLB = fn_Tind(conc{2},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{2},TpUB,ProbeLocXY(1),WaveType);
        
        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time + T_windowSec;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time - T_windowSec ;
        
        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate); %fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn2 = mean(Sn_mat(:,jj0:jj1), 2);

        m0Array = trapz(ff, Sn_mat);
        Hs2(:,i) = 4*sqrt(m0Array);
        Tm2(:,i) = tm_vec - t0;    
        TmE2(1,i) = 0;
        TmE2(2,i) = t1 - t0;

        ff2{i} = ff;

        SnS2{i} = Sn2;
        Times2(i,:) = [tm_vec(jj0) - t0 ;tm_vec(jj1) - t0];
        Times2US(i,:) = [tm_vec(jj0);tm_vec(jj1)];
        i = i +1;
    end

    SnS1L = [];
    SnS2L = [];
    for i = 1: length( SnS1)
        SnS1L = [SnS1L , SnS1{i} ];
        SnS2L = [SnS2L , SnS2{i} ];
    end

    AverageTimeWindows1 = mean(SnS1L,2);
    AverageTimeWindows2 = mean(SnS2L,2);

    %Jonswap
    wJS = 2*pi/Tp;
    ffJS = 0:0.01:5;
    JSspec  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
   
    figure();
    hold on;
    Description = ['Concentations (', num2str(conc{1}) ,' , ', num2str(conc{2}) ') Tp = ', num2str(Tp), '[s] Hs = ', num2str(Hs) , '[m] Fourier Interval ', num2str(FourierWindow_Sec),'[s] Time Interval : ( ',num2str(Times2US(1,1)),'[s] , ', num2str(Times2US(1,2)),'[s] )' ];
    for i = 1: length( SnS1)
        if i == 1
            p1l = plot(ff1{i}, SnS1{i}, '-r', 'DisplayName', ['Conc ', num2str(conc{1})  ' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )']);
            hold on;
            p2l = plot(ff2{i}, SnS2{i}, '-b', 'DisplayName', ['Conc ', num2str(conc{2})  ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )']);
        end
            plot(ff1{i}, SnS1{i}, '-r')
            plot(ff2{i}, SnS2{i}, '-b')
    end
    p3 = plot(ffJS,JSspec, '--k', 'DisplayName', 'JonSwap Tp, Hs', 'LineWidth', 2);   
    xlabel('Frequency [Hz]');
    ylabel('Spectra');
    title(['Individual Probes (Average over Time Intervals) ', Description]);
    legend([p1l,p2l, p3]);
    xlim([0 1.0/ (Tp/4)])
    hold off

    
    %Average over Probes Reflection vs Transmission
    figure();
    plot(ff1{1}, AverageTimeWindows1, '-r', 'DisplayName',  ['Conc ', num2str(conc{1})  ' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(ff2{1}, AverageTimeWindows2, '-b', 'DisplayName',  ['Conc ', num2str(conc{2})  ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(ffJS,JSspec, '--k', 'DisplayName', 'JonSwap Tp, Hs', 'LineWidth', 2);   
    xlabel('Frequency [Hz]');
    ylabel('Spectra');
    title([' Average Spectra (Probes and Time Windows) ', Description]);
    legend();
    xlim([0 1.0/ (Tp/4)])
    hold off;
    
    figure();
    for i = 1: size(TmE1,2)
        
        if i == 1
           p1 = plot(Tm1(:,i), Hs1(:,i), '-b', 'DisplayName','Reflection Probes', 'LineWidth', 2);
           hold on;
           p2 = plot(Tm2(:,i), Hs2(:,i), '-r', 'DisplayName','Transmission Probes', 'LineWidth', 2);
        end
           plot(Tm1(:,i), Hs1(:,i), '-b', 'LineWidth', 2);
           plot(Tm2(:,i), Hs2(:,i), '-r', 'LineWidth', 2);        
    end
    
    MinTimes1 = max(max(min(Times1)),max(min(Times2)));
    MaxTimes1 = min(min(max(Times1)),min(max(Times2)));
    
    p3 = plot([min(min(Tm1)),max(max(Tm1))], [Hs, Hs], '--k', 'DisplayName','Hs', 'LineWidth', 2);
    p4 = plot([MinTimes1,MinTimes1], [0, 2*Hs], '-g', 'DisplayName', 'Time Window Interval', 'LineWidth', 2);
    plot([MaxTimes1,MaxTimes1], [0, 2*Hs], '-g', 'LineWidth', 2);
    hold off;
    title(['Significant Wave Heights ', Description ])
    xlabel('Average Time (Translated To First Waves Reach) (s)');
    ylabel('Significant Wave Height (m)');
    legend([p1 p2 p3 p4]);

end



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

