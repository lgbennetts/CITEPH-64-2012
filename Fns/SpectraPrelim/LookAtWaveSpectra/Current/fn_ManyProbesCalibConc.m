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

function fn_ManyProbesCalibConc(conc,TestName,Tp,probes)

%Compare calibration to disk
% if ~exist('conc','var') conc={0,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'16','8'}}; end ;  % 9, 14
% {{'15','10'}, {'14','9'},{'16','8'}};

if ~exist('conc','var') conc={0,79}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'17','22'}}; end ;  % 9, 14
%{{'14','19'}, {'15','20'},{'16','21'},{'17','22'}};

%'17' '22'

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
figure();
hold on;
for j = 1: length(TestName)
    

    t0description ='waves reach x';
    t1description ='final waves reach x';

    for i = 1:length(probes{1})
     %Read Data in Calib 1, at probe
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, probes{1}(i));

        Tp = c_pram.period/10.0;
        dt = tm(2) - tm(1);
        
        eval(TpersStr);
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

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        %% Fourier Transform

        ff1{i} = ff;

        SnS1{i} = Sn1;
        Times1{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];



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
        
        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn2 = mean(Sn_mat(:,jj0:jj1), 2);

        %% Fourier Transform
%         [Sn2,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
%         period = 1./ff;


        ff2{i} = ff;

        SnS2{i} = Sn2;
        Times2{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];

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
    hmJS = Hs;
%     hmJS = Hs/10.0;
    ffJS = 0:0.01:5;
    JSspec  = jonswap(2*pi.*ffJS,'wp',wJS,'Hs',hmJS);
    JSspec1  = jonswap(2*pi.*ff1{1},'wp',wJS,'Hs',hmJS);
    
    %get transmission coefficients as function of k.
    
    %Models
    if conc{2}==39
    dum_c = 100*pi*(0.495^2)/2;
    elseif conc{2}==79
    dum_c = 100*pi*(0.495^2);
    end
    
    %Predictions
    Boltzman0 = Main_AttnModels(model_pers,dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
    Boltzman0 = (Boltzman0.value).^2;
    Boltzman1 = Main_AttnModels(model_pers,dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
    Boltzman1 = (Boltzman1.value).^2;
    %2Demm prediction
    TwoDEMM = Main_AttnModels(model_pers,dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
    TwoDEMM = (TwoDEMM.value).^2;

    PTC_2EEM = interp1(1./model_pers,TwoDEMM,ff1{1});
    PTC_Bolt0 = interp1(1./model_pers,Boltzman0,ff1{1});
    PTC_Bolt1 = interp1(1./model_pers,Boltzman1,ff1{1});

    
    PTC_2EEMJS = interp1(1./model_pers,TwoDEMM,ffJS);
    PTC_Bolt0JS = interp1(1./model_pers,Boltzman0,ffJS);
    PTC_Bolt1JS = interp1(1./model_pers,Boltzman1,ffJS);
    
    Description = ['Concentations (', num2str(conc{1}) ,' , ', num2str(conc{2}) ') Tp = ', num2str(Tp), '[s] Hs = ', num2str(Hs) , '[m] Fourier Interval ', num2str(FourierWindow_Sec),'[s] Time Interval : ( ',num2str(Times2{1}(1)),'[s] , ', num2str(Times2{1}(2)),'[s] )' ];
    for i = 1: length( SnS1)
        if i == 1
            p1l = plot(ff1{i}, SnS1{i}, '-r', 'DisplayName', ['Conc ', num2str(conc{1})  ' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )']);
            hold on;
            p2l = plot(ff2{i}, SnS2{i}, '-b', 'DisplayName', ['Conc ', num2str(conc{2})  ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )']);
        end
            plot(ff1{i}, SnS1{i}, '-r')
            plot(ff2{i}, SnS2{i}, '-b')
    end
    xlabel('Frequency [Hz]');
    ylabel('Spectra');
    title(['Individual Probes (Average over Time Intervals) ', Description]);
    legend([p1l,p2l]);
    xlim([0 1.0/ (Tp/4)])
    hold off
    
    %Plot of spectra and predictions
    
    %Compare Calib and Conc
    figure();
    plot(ff1{1}, AverageTimeWindows1, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(ffJS, JSspec, '-k', 'DisplayName',  'JonSwap' , 'LineWidth', 2);
    plot(ff1{1}, AverageTimeWindows2, '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(ff1{1}, (PTC_2EEM).*AverageTimeWindows1, '--k', 'DisplayName', '2DEMM Predicted Spectra' , 'LineWidth', 2);
    plot(ff1{1}, (PTC_Bolt0).*AverageTimeWindows1, '-.r', 'DisplayName', 'Boltzman Predicted Spectra' , 'LineWidth', 2);
    plot(ff1{1}, (PTC_Bolt1).*AverageTimeWindows1, '-.g', 'DisplayName', 'Boltzman (Remove Scattering) Predicted Spectra' , 'LineWidth', 2);
    xlabel('Frequency [Hz]');
    ylabel('Spectra ');
    title(['Average Spectra (Probes and Time Intervals) with Predicted Spectra ', Description]);
    legend();
    xlim([0 1.0/ (Tp/4)]);
    hold off;

    figure();
    plot(ffJS, JSspec, '-r', 'DisplayName',  'JonSwap' , 'LineWidth', 2);
    hold on;
    plot(ff1{1}, AverageTimeWindows2, '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(ffJS, PTC_2EEMJS.*JSspec, '--k', 'DisplayName', '2DEMM Predicted Spectra' , 'LineWidth', 2);
    plot(ffJS, PTC_Bolt0JS.*JSspec, '-.r', 'DisplayName', 'Boltzman Predicted Spectra' , 'LineWidth', 2);
    plot(ffJS, PTC_Bolt1JS.*JSspec, '-.g', 'DisplayName', 'Boltzman (Remove Scattering) Predicted Spectra' , 'LineWidth', 2);
    xlabel('Frequency [Hz]');
    ylabel('Spectra ');
    title(['Average Spectra (Probes and Time Intervals) with Predicted JONSWAP Spectra ', Description]);
    legend();
    xlim([0 1.0/ (Tp/4)]);
    hold off;
    
    %Transmission Coefficients
    figure()
    plot(1./ff1{1}, AverageTimeWindows2 ./ AverageTimeWindows1,'xr','DisplayName','Data' );
    hold on;
    plot(model_pers,TwoDEMM,'--k','DisplayName','2D EEM' )
    plot(model_pers,Boltzman0,'--b','DisplayName','Boltzman' )
    plot(model_pers,Boltzman1,'--r','DisplayName','Boltzman (Remove Scattering)' )
    xlabel('Periods');
    ylabel('Transmission Coefficient');
    title(['Transmission Coefficient', Description]);
    legend();
    xlim([0,3]);  
    ylim([0,1]);
    hold off;
    
    figure()
    plot(1./ff1{1}, AverageTimeWindows2 ./ JSspec1,'xr','DisplayName','Data' );
    hold on;
    plot(model_pers,TwoDEMM,'--k','DisplayName','2D EEM' )
    plot(model_pers,Boltzman0,'--b','DisplayName','Boltzman' )
    plot(model_pers,Boltzman1,'--r','DisplayName','Boltzman (Remove Scattering)' )
    xlabel('Periods (s)');
    ylabel('Transmission Coefficient');
    title(['Transmission Coefficient JonSwap', Description]);
    legend();
    xlim([0,3]);  
    ylim([0,1]);
    hold off;
    
 
    

end



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

