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
% if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  

if ~exist('conc','var') conc={1,79}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'14','19'}}; end ;  % {'14','9'}, {'15','10'}
%{'14','19'}

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
close all;
t0description ='waves reach x';
t1description ='final waves reach x';

for j = 1: length(TestName)
 
    
    %Tp
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{1}{1}, probes{1}(1));

    Tp = c_pram.period/10.0;
%     TpFac = [3/4,10/12,11/12,12/12,13/12,14/12,15/12,16/12,17/12,18/12];
    TpFac = [1/12,3/12,5/12,7/12,9/12,11/12,13/12,15/12,17/12,19/12,21/12,23/12,25/12];

    S1A = [];
    S2A = [];
    PA = [];
    
    for k = 1: length(TpFac)-1
        LBfactor = TpFac(k);
        UBfactor = TpFac(k+1);

        CFactor = 0.5*(LBfactor +UBfactor);


        for i = 1:length(probes{1})
         %Read Data in Calib 1, at probe
            [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, probes{1}(i));

            TpC = CFactor*Tp;
            TpLB = Tp*LBfactor;
            TpUB = Tp*UBfactor;
            
            Tpers = 100;
            T_window = Tpers*TpC;

            Hs = c_pram.wave_height/100.0;

            dt = tm(2) - tm(1);

            Tind = fn_Tind(conc{1},TpC,ProbeLocXY(1),WaveType); %Slower
            TindLB = fn_Tind(conc{1},TpLB,ProbeLocXY(1),WaveType);
            TindUB = fn_Tind(conc{1},TpUB,ProbeLocXY(1),WaveType);


            t0index = find(strcmp({TindLB.description},t0description));
            t0 = TindLB(t0index).time + T_window;
            t1index = find(strcmp({TindUB.description},t1description));
            t1 = TindUB(t1index).time - T_window ;

            [Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,TpC,T_window);
            period = 1./ff;

            [~,jj0]=min(abs(tm_vec-t0));
            [~,jj1]=min(abs(tm_vec-t1));

            Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
            %% Fourier Transform

            Periods1{i} = period;

            SnS1{i} = Sn1;
            Times1{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];



            [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{j}{2}, probes{2}(i));

            dt = tm(2) - tm(1);

            Tind = fn_Tind(conc{2},TpC,ProbeLocXY(1),WaveType); %Slower
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
        
        %Get Periods between limit and spectra
        
        S1C = AverageTimeWindows1(and(Periods1{1} > TpLB, Periods1{1} < TpUB ));
        S2C = AverageTimeWindows2(and(Periods1{1} > TpLB, Periods1{1} < TpUB ));
        PC =  Periods1{1}(and(Periods1{1} > TpLB, Periods1{1} < TpUB ));
        
        S1A = cat(1,S1A,flip(S1C));
        S2A = cat(1,S2A,flip(S2C));
        PA = cat(1,PA,flip(PC));
    end

    %Area under curve
    Cff = 1./PA;
    CMS = S1A;
    m0 = abs(trapz(Cff, CMS));
    
    
    %Jonswap
    wJS = 2*pi/Tp;
    hmJS = 4*sqrt(m0);
%     hmJS = Hs/10.0;
    JSspec  = abs(jonswap(2*pi./PA,'wp',wJS,'Hs',hmJS));

    
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
    
    PredictedTransmissionCoefficient = interp1(model_pers,TwoDEMM,PA);

    Description = ['Concentations (', num2str(conc{1}) ,' , ', num2str(conc{2}) ') Tp = ', num2str(Tp), '[s] Hs = ', num2str(Hs) , '[m] Fourier Interval ', num2str(T_window),'[s] Time Interval : ( ',num2str(Times2{1}(1)),'[s] , ', num2str(Times2{1}(2)),'[s] )' ];
    

    %No Moving Average
    figure();
    plot(PA, S1A, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(PA, S2A, '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(PA,JSspec,'--k', 'LineWidth', 2, 'DisplayName', ['JonSwap wp: ', num2str(wJS), '[Hz]' , '  Hs (4sqrt(m0)) ', num2str(round(hmJS,4)), '[m]' ])
    xlabel('Periods');
    ylabel('Spectra');
    title(['Average Spectra (Probes and Time Intervals) ', Description]);
    legend();
    xlim([0,2*Tp]);
    hold off
    
    %Moving Average
    m = 7;
    figure();
    plot(PA, movmean(S1A,m), '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(PA, movmean(S2A,m), '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(PA,JSspec,'--k', 'LineWidth', 2, 'DisplayName', ['JonSwap wp: ', num2str(wJS), '[Hz]' , '  Hs (4sqrt(m0)) ', num2str(round(hmJS,4)), '[m]' ])
    xlabel('Periods');
    ylabel('Spectra');
    title(['Moving Average (',num2str(m),')Average Spectra (Probes and Time Intervals) ', Description]);
    legend();
    xlim([0,2*Tp]);
    hold off
    
    figure();
    plot(PA, S1A, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(PA, S2A, '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(PA, PredictedTransmissionCoefficient.*S1A, '-k', 'DisplayName', '2DEMM Predicted Spectra (From Calibration Spectra)' , 'LineWidth', 2);
    xlabel('Periods');
    ylabel('Spectra');
    title(['Average Spectra (Probes and Time Intervals) ', Description]);
    legend();
    xlim([0,2*Tp]);
    hold off
    %Moving Average
    

    figure();
    plot(PA, movmean(S1A,m), '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )'] , 'LineWidth', 2);
    hold on;
    plot(PA, movmean(S2A,m), '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
    plot(PA,PredictedTransmissionCoefficient.*movmean(S1A,m),'-k', 'LineWidth', 2, 'DisplayName', '2DEMM Predicted Spectra (From Calibration Spectra)')
    xlabel('Periods');
    ylabel('Spectra');
    title(['Moving Average (',num2str(m),')Average Spectra (Probes and Time Intervals) ', Description]);
    legend();
    xlim([0,2*Tp]);
    hold off
    
    
    %Figure - Transmission coefficent
    
    %Models
    if conc{2}==39
    dum_c = 100*pi*(0.495^2)/2;
    elseif conc{2}==79
    dum_c = 100*pi*(0.495^2);
    end

    Boltzman0 = Main_AttnModels(model_pers,dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
    Boltzman0 = (Boltzman0.value).^2;
    Boltzman1 = Main_AttnModels(model_pers,dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
    Boltzman1 = (Boltzman1.value).^2;

    TwoDEMM = Main_AttnModels(model_pers,dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
    TwoDEMM = (TwoDEMM.value).^2;

    figure();
    plot(model_pers,Boltzman0,'--r','DisplayName','Boltzmann Remove Scattering','LineWidth',2);
    hold on;
    plot(model_pers,Boltzman1,'--b','DisplayName','Boltzmann','LineWidth',2);

    plot(model_pers,TwoDEMM,'-.k','DisplayName','2D EMM','LineWidth',2);

    TransCoeff = S2A ./ S1A; 
    plot(PA,TransCoeff, '.k' );
    
    xlabel('Periods');
    ylabel('Transmission Coefficient');
    title(['Transmission Coefficient', Description]);    
    xlim([0,2*Tp]);
    ylim([0,1]);
    hold off;

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

