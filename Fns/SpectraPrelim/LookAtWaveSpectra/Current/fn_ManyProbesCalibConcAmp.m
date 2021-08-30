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

function fn_CompareSpectraCalibConcOverlayAverageMultiple(conc,TestName,Tp,probes)

%Compare calibration to disk
% if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  

if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  % 9, 14


if ~exist('Cols','var'); Cols= {'xr', '^b','sk'}; end ;

if ~exist('wbinwidth','var');  wbinwidth=6; end;
if ~exist('UBfactor','var');  UBfactor=1.5; end;
if ~exist('LBfactor','var');  LBfactor=0.75; end;

if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
% close all;

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
        
        Tpers = 100;
        T_window = Tpers*Tp;
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs = c_pram.wave_height/100.0;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        TindLB = fn_Tind(conc{1},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{1},TpUB,ProbeLocXY(1),WaveType);


        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time + T_window;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time - T_window ;

        [Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        %% Fourier Transform

        Periods1{i} = period;

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
        t0 = TindLB(t0index).time + T_window;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time - T_window ;
        
        [Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn2 = mean(Sn_mat(:,jj0:jj1), 2);

        %% Fourier Transform
%         [Sn2,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
%         period = 1./ff;


        Periods2{i} = period;

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

    %Area under curve
    Cff = 1./Periods1{i}(and(Periods1{1} > Tp/3, Periods1{1} < 3*Tp));
    CMS = AverageTimeWindows1(and(Periods1{1} > Tp/3, Periods1{1} < 3*Tp));
    m0 = trapz(Cff, CMS);
    
    %Jonswap
    wJS = 2*pi/Tp;
    hmJS = 4*sqrt(m0);
%     hmJS = Hs/10.0;
    JSspec  = jonswap(2*pi./Periods1{1},'wp',wJS,'Hs',hmJS);

    
    %get transmission coefficients as function of k.
    
    %Models
    if conc{2}==39
    dum_c = 100*pi*(0.495^2)/2;
    elseif conc{2}==79
    dum_c = 100*pi*(0.495^2);
    end
    
    %2Demm prediction
    TwoDEMM = Main_AttnModels(model_pers,dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
    TwoDEMM = (TwoDEMM.value).^2;
    
    PredictedTransmissionCoefficient = interp1(model_pers,TwoDEMM,Periods1{1});


    Description = ['Concentations (', num2str(conc{1}) ,' , ', num2str(conc{2}) ') Tp = ', num2str(Tp), '[s] Hs = ', num2str(Hs) , '[m] Fourier Interval ', num2str(T_window),'[s] Time Interval : ( ',num2str(Times2{1}(1)),'[s] , ', num2str(Times2{1}(2)),'[s] )' ];
    for i = 1: length( SnS1)
        if i == 1
            p1l = plot(Periods1{i}, SnS1{i}, '-r', 'DisplayName', ['Conc ', num2str(conc{1})  ' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )']);
            hold on;
            p2l = plot(Periods2{i}, SnS2{i}, '-b', 'DisplayName', ['Conc ', num2str(conc{2})  ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )']);
        end
            plot(Periods1{i}, SnS1{i}, '-r')
            plot(Periods2{i}, SnS2{i}, '-b')
    end
    xlabel('Periods');
    ylabel('Spectra');
    title(['Individual Probes (Average over Time Intervals) ', Description]);
    legend([p1l,p2l]);
    xlim([0,2*Tp]);
    hold off
    
    figure();
    plot(Periods1{i}, AverageTimeWindows1, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(Periods2{i}, AverageTimeWindows2, '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
%     plot(Periods1{i},JSspec,'--k', 'LineWidth', 2, 'DisplayName', ['JonSwap wp: ', num2str(wJS), '[Hz]' , '  Hs (4sqrt(m0)) ', num2str(round(hmJS,4)), '[m]' ])
    plot(Periods1{1}, PredictedTransmissionCoefficient.*AverageTimeWindows1, '--k', 'DisplayName', '2DEMM Predicted Spectra' , 'LineWidth', 2);
    

    xlabel('Periods');
    ylabel('Spectra (Amplitudes sqrt(2*S))');
    title(['Average Spectra (Probes and Time Intervals) ', Description]);
    legend();
    xlim([0,2*Tp]);
    hold off
    
    %Moving average - 7
    m = 7;
    MovAvg1 =  movmean(AverageTimeWindows1,m);
    MovAvg2 =  movmean(AverageTimeWindows2,m);
    figure();
    plot(Periods1{i}, MovAvg1, '-r', 'DisplayName',  ['Conc ', num2str(conc{1})] , 'LineWidth', 2);
    hold on;
    plot(Periods2{i}, MovAvg2, '-b', 'DisplayName',  ['Conc ', num2str(conc{2})] , 'LineWidth', 2)    ;
    plot(Periods1{1}, PredictedTransmissionCoefficient.*MovAvg1, '--k', 'DisplayName', '2DEMM Predicted Spectra' , 'LineWidth', 2);
    
    xlabel('Periods');
    ylabel('Spectra');
    title(['Moving Average (',num2str(m),') of Average Spectra (Probes and Time Intervals) ', Description]);
    legend();
    xlim([0,2*Tp]);
    hold off;
%     
%     %Figure - Transmission coefficent
%     
%     %Models
%     if conc{2}==39
%     dum_c = 100*pi*(0.495^2)/2;
%     elseif conc{2}==79
%     dum_c = 100*pi*(0.495^2);
%     end

%     Boltzman0 = Main_AttnModels(model_pers,dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
%     Boltzman0 = (Boltzman0.value).^2;
%     Boltzman1 = Main_AttnModels(model_pers,dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
%     Boltzman1 = (Boltzman1.value).^2;
% 
%     TwoDEMM = Main_AttnModels(model_pers,dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
%     TwoDEMM = (TwoDEMM.value).^2;
% 
%     figure();
%     plot(model_pers,Boltzman0,'--r','DisplayName','Boltzmann Remove Scattering','LineWidth',2);
%     hold on;
%     plot(model_pers,Boltzman1,'--b','DisplayName','Boltzmann','LineWidth',2);
% 
%     plot(model_pers,TwoDEMM,'-.k','DisplayName','2D EMM','LineWidth',2);

%     TransCoeff = AverageTimeWindows2 ./ AverageTimeWindows1; 
%     plot(Periods1{i}(and (Periods1{i} > TpLB ,Periods1{i} < TpUB )),TransCoeff(and (Periods1{i} > TpLB ,Periods1{i} < TpUB )), '.k' );
%     
%     TransCoeffMov = movsum(AverageTimeWindows2,wbinwidth) ./ movsum(AverageTimeWindows1,wbinwidth); 
%     plot(Periods1{i}(and (Periods1{i} > TpLB ,Periods1{i} < TpUB )),TransCoeffMov(and (Periods1{i} > TpLB ,Periods1{i} < TpUB )), '.r' );
% 
%     xlabel('Periods');
%     ylabel('Transmission Coefficient');
%     title(['Transmission Coefficient', Description]);    
%     xlim([0,2*Tp]);
%     ylim([0,1]);
%     hold off;

    %
    

end



return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

