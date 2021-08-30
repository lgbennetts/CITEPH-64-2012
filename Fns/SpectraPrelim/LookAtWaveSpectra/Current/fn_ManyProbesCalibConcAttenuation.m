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
% if ~exist('conc','var') conc={0,39,0,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20,11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'16','8','15','10'}}; end ;  % 9, 14


if ~exist('conc','var') conc={0,79,0,79}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20,11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'16','21','15','20'}}; end ;  % 9, 14

%{{'16','8','15','10'}}

% if ~exist('conc','var') conc={0,79}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'17','22'}}; end ;  % 9, 14
% %{{'14','19'}, {'15','20'},{'16','21'},{'17','22'}};

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
if ~exist('model_pers','var'); model_pers=[10^-12,0.3:0.02:2,10^3]; end
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

        Tp = c_pram.period/10.0;
        dt = tm(2) - tm(1);
        
        eval(TpersStr);
        T_windowSec = round(Tpers*Tp);
        SamplingRate = floor(1./dt);
        T_window = SamplingRate*T_windowSec;
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs12 = c_pram.wave_height/100.0;


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
        ff2{i} = ff;

        SnS2{i} = Sn2;
        Times2{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];
        
        
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{3},TestName{j}{3}, probes{3}(i));

        Tp = c_pram.period/10.0;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{3},Tp,ProbeLocXY(1),WaveType); %Slower
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs34 = c_pram.wave_height/100.0;
        
        TindLB = fn_Tind(conc{3},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{3},TpUB,ProbeLocXY(1),WaveType);
        
        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time + T_windowSec;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time - T_windowSec ;
        
        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn3 = mean(Sn_mat(:,jj0:jj1), 2);
        ff3{i} = ff;

        SnS3{i} = Sn3;
        Times3{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];

        
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{4},TestName{j}{4}, probes{4}(i));

        Tp = c_pram.period/10.0;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{4},Tp,ProbeLocXY(1),WaveType); %Slower
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        
        TindLB = fn_Tind(conc{4},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{4},TpUB,ProbeLocXY(1),WaveType);
        
        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time + T_windowSec;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time - T_windowSec ;
        
        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn4 = mean(Sn_mat(:,jj0:jj1), 2);
        ff3{i} = ff;

        SnS4{i} = Sn4;
        Times4{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];
        i = i +1;
    end

    SnS1L = [];
    SnS2L = [];
    SnS3L = [];
    SnS4L = [];
    for i = 1: length( SnS1)
        SnS1L = [SnS1L , SnS1{i} ];
        SnS2L = [SnS2L , SnS2{i} ];
        SnS3L = [SnS3L , SnS3{i} ];
        SnS4L = [SnS4L , SnS4{i} ];
    end

    ATW1 = mean(SnS1L,2);
    ATW2 = mean(SnS2L,2);
    ATW3 = mean(SnS3L,2);
    ATW4 = mean(SnS4L,2);
    
    ATW1max = max(SnS1L,[],2);
    ATW1min = min(SnS1L,[],2);
    ATW2max = max(SnS2L,[],2);
    ATW2min = min(SnS2L,[],2);
    ATW3max = max(SnS3L,[],2);
    ATW3min = min(SnS3L,[],2);
    ATW4max = max(SnS4L,[],2);
    ATW4min = min(SnS4L,[],2);
    
    
    %Jonswap
    wJS = 2*pi/Tp;
    hmJS = Hs12;
%     hmJS = Hs/10.0;
    ffJS = min(ff1{1}):0.01:max(ff1{1});
    JSspec12  = jonswap(2*pi.*ffJS,'wp',wJS,'Hs',Hs12);
    JSspec34  = jonswap(2*pi.*ffJS,'wp',wJS,'Hs',Hs34);
    
    %get transmission coefficients as function of k.
    
    %Models
    if conc{2}==39
    dum_c = 100*pi*(0.495^2)/2;
    elseif conc{2}==79
    dum_c = 100*pi*(0.495^2);
    end
    
    %Predictions
    Boltzman0 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
    Boltzman0 = (Boltzman0.value).^2;
    Boltzman0 = Boltzman0.';
    Boltzman0 = [0;Boltzman0;1];
    
    
    Boltzman1 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
    Boltzman1 = (Boltzman1.value).^2;
    Boltzman1 = Boltzman1.';
    Boltzman1 = [0;Boltzman1;1];
    
    %2Demm prediction
    TwoDEMM = Main_AttnModels(model_pers(2:end-1),dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
    TwoDEMM = (TwoDEMM.value).^2;
    TwoDEMM = TwoDEMM.';
    TwoDEMM = [0;TwoDEMM;1];
    
    ffp = ff1{1}*Tp;
    ffJSp = (ffJS*Tp);

    
    Description = ['Concentations (', num2str(conc{1}) ,' , ', num2str(conc{2}) ') Tp = ', num2str(Tp), '[s] Fourier Interval ', num2str(FourierWindow_Sec),'[s] Time Interval : ( ',num2str(Times2{1}(1)),'[s] , ', num2str(Times2{1}(2)),'[s] )' ];

    %Compare calibration tests
    m01 = trapz(ff1{1},ATW1);
    m02 = trapz(ff1{1},ATW2);
    m03 = trapz(ff1{1},ATW3);
    m04 = trapz(ff1{1},ATW4);
    
    HS1 = 4*sqrt(m01);
    HS2 = 4*sqrt(m02);
    HS3 = 4*sqrt(m03);
    HS4 = 4*sqrt(m04);

    
    figure();
    plot(ffp, ATW1./m01, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 2);
    hold on;
    plot(ffp, ATW3./m03, '-b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Probes : ( ', num2str(probes{3}(1)),' , ', num2str(probes{3}(end)) ' ) | Hs ', num2str(Hs34),' [m]'] , 'LineWidth', 2);
    plot(ffJSp, JSspec12./(Hs12/4)^2, '--k', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 2);
    xlabel('frequency / fp');
    ylabel('Spectra / m0 ');
    title(['Calibration Spectra Comparisons', Description]);
    legend();
    set(gca, 'YScale', 'log')
    ylim([10^-12 , 10 ]);
    hold off; 

    
    %Figure
    figure();
    p1 = plot(ffp, ATW1./m01, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 3);
    hold on;
    ATWmaxp = ATW1max./m01;
    ATWminp = ATW1min./m01;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    alpha(0.3)
    p4 = plot(ffp, ATWmaxp, '--r', 'LineWidth', 1, 'DisplayName','max and min');
    plot(ffp, ATWminp, '--r', 'LineWidth', 1);
    
    p3 = plot(ffp, ATW3./m03, '-b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Probes : ( ', num2str(probes{3}(1)),' , ', num2str(probes{3}(end)) ' ) | Hs ', num2str(Hs34),' [m]'] , 'LineWidth', 3);
    ATWmaxp = ATW3max./m03;
    ATWminp = ATW3min./m03;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    alpha(0.3)
    p5 = plot(ffp, ATWmaxp, '--b', 'LineWidth', 1, 'DisplayName','max and min');
    plot(ffp, ATWminp, '--b', 'LineWidth', 1);   
    p2 = plot(ffJSp, JSspec12./(Hs12/4)^2, '--k', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 3);
    xlabel('frequency / fp');
    ylabel('Spectra / m0 ');
    title(['Calibration Spectra Comparisons With Error ', Description]);
    legend([p1 p2 p3 p4 p5]);
    set(gca, 'YScale', 'log')
    ylim([10^-12 , 10 ]);
    hold off; 
    
    %Figure - Spectra Comparisons
    figure();
    p2 = plot(ffJSp, JSspec12./(Hs12/4)^2, '-k', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 3);
    hold on;
    p1 = plot(ffp, ATW1./m01, '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 3);
    ATWmaxp = ATW1max./m01;
    ATWminp = ATW1min./m01;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    alpha(0.3)
    p4 = plot(ffp, ATWmaxp, ':r', 'LineWidth', 1, 'DisplayName','+- 2 std');
    plot(ffp, ATWminp, ':r', 'LineWidth', 1);
    p3 = plot(ffp, ATW3./m03, '-b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Probes : ( ', num2str(probes{3}(1)),' , ', num2str(probes{3}(end)) ' ) | Hs ', num2str(Hs34),' [m]'] , 'LineWidth', 3);
    ATWmaxp = ATW3max./m03;
    ATWminp = ATW3min./m03;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    alpha(0.3)
    p5 = plot(ffp, ATWmaxp, ':b', 'LineWidth', 1, 'DisplayName','+- 2 std');
    plot(ffp, ATWminp, ':b', 'LineWidth', 1); 

    p6 = plot(ffp, ATW2./m02, '--r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 3);
    ATWmaxp = ATW2max./m02;
    ATWminp = ATW2min./m02;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    alpha(0.3)
    p7 = plot(ffp, ATWmaxp, '-.r', 'LineWidth', 1, 'DisplayName','+- 2 std');
    plot(ffp,ATWminp, '-.r', 'LineWidth', 1);

    p8 = plot(ffp, ATW4./m04, '--b', 'DisplayName',  ['Conc ', num2str(conc{4}) ,' | Probes : ( ', num2str(probes{4}(1)),' , ', num2str(probes{4}(end)) ' ) | Hs ', num2str(Hs34),' [m]'] , 'LineWidth', 3);
    ATWmaxp = ATW4max./m04;
    ATWminp = ATW4min./m04;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    alpha(0.3)
    p9 = plot(ffp, ATWmaxp, '-.b', 'LineWidth', 1, 'DisplayName','+- 2 std');
    plot(ffp, ATWminp, '-.b', 'LineWidth', 1);
    
    xlabel('frequency / fp');
    ylabel('Spectra / m0 ');
    title(['Calibration Spectra Comparisons With Error ', Description]);
    legend([p1 p2 p3 p4 p5 p6 p7 p8 p9]);
    set(gca, 'YScale', 'log')
    ylim([10^-12 , 10 ]);
    hold off; 
    
    
    ffmod = 1.0./ model_pers.';
    ffmodp  = ffmod(2:end-1).*Tp;
    
    Amp1 = sqrt(2.*ATW1);
    Amp3 = sqrt(2.*ATW3);
    
    Tran21 = (ATW2./ATW1);
    Tran21Max = (ATW2max) ./ (ATW1);
    Tran21Min = (ATW2min) ./ (ATW1);
    Tran43 = (ATW4./ATW3);
    Tran43Max = (ATW4max) ./ (ATW3);
    Tran43Min = (ATW4min) ./ (ATW3);
    
    TEM = TwoDEMM(2:end-1);
    TB0 = Boltzman0(2:end-1);
    TB1 = Boltzman1(2:end-1);
    
    %Figure
    figure();
    p1 = plot(ffp, Tran21, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
    hold on;
    ATWmaxp = Tran21Max;
    ATWminp = Tran21Min;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    alpha(0.3)
    p1a = plot(ffp, ATWmaxp, ':r', 'LineWidth', 1, 'DisplayName','max/min');
    plot(ffp, ATWminp, ':r', 'LineWidth', 1); 
    
    p2 = plot(ffp, Tran43, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
    ATWmaxp = Tran43Max;
    ATWminp = Tran43Min;
    ShadeX = [ffp; flipud(ffp)];
    inBetween = [(ATWmaxp); flipud(ATWminp)];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    alpha(0.3)
    p2a = plot(ffp, ATWmaxp, ':b', 'LineWidth', 1, 'DisplayName','max/min');
    plot(ffp, ATWminp, ':b', 'LineWidth', 1);     
    
    
    p3 = plot(ffmodp,TEM,'--k', 'DisplayName','2DEMM', 'LineWidth',3);
    p4 = plot(ffmodp,TB0,'-.b', 'DisplayName','Boltzmann', 'LineWidth',3);
    p5 = plot(ffmodp,TB1,'-.r', 'DisplayName','Boltzmann Remove Scattering', 'LineWidth',3);
    xlabel('Frequency / Fp');
    ylabel('Transmission Coefficient');  
     title(['Transmission Coefficient With Error ', Description]);
    xlim([0,4]);
    legend([p1 p1a p2 p2a p3 p4 p5]);

    %limit to 1
%     Tran21(Tran21 > 1) = 1;
%     Tran21Max(Tran21Max > 1) = 1;
%     Tran21Min(Tran21Min > 1) = 1;
%     Tran43(Tran43 > 1) = 1;
%     Tran43Max(Tran43Max > 1) = 1;
%     Tran43Min(Tran43Min > 1) = 1;   
%     figure();
%     p1 = plot(ffp, Tran21, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
%     hold on;
%     ATWmaxp = Tran21Max;
%     ATWminp = Tran21Min;
%     ShadeX = [ffp; flipud(ffp)];
%     inBetween = [(ATWmaxp); flipud(ATWminp)];
%     patch(ShadeX, inBetween, 'r', 'LineStyle','none')
%     alpha(0.3)
%     p1a = plot(ffp, ATWmaxp, ':r', 'LineWidth', 1, 'DisplayName','max/min');
%     plot(ffp, ATWminp, ':r', 'LineWidth', 1); 
%     
%     p2 = plot(ffp, Tran43, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
%     ATWmaxp = Tran43Max;
%     ATWminp = Tran43Min;
%     ShadeX = [ffp; flipud(ffp)];
%     inBetween = [(ATWmaxp); flipud(ATWminp)];
%     patch(ShadeX, inBetween, 'b', 'LineStyle','none')
%     alpha(0.3)
%     p2a = plot(ffp, ATWmaxp, ':b', 'LineWidth', 1, 'DisplayName','max/min');
%     plot(ffp, ATWminp, ':b', 'LineWidth', 1);     
%     
%     
%     p3 = plot(ffmodp,TEM,'--k', 'DisplayName','2DEMM', 'LineWidth',3);
%     p4 = plot(ffmodp,TB0,'-.b', 'DisplayName','Boltzmann', 'LineWidth',3);
%     p5 = plot(ffmodp,TB1,'-.r', 'DisplayName','Boltzmann Remove Scattering', 'LineWidth',3);
%     xlabel('Frequency / Fp');
%     ylabel('Transmission Coefficient');  
%     title(['Transmission Coefficient Max 1 ', Description]);
%     xlim([0,4]);
%     legend([p1 p1a p2 p2a p3 p4 p5]);
    
    
    %attenuation
    AEM = -(log(TEM)./5);
    AB1 = -(log(TB1)./5);
    AB0 = -(log(TB0)./5);
    A21 = -(log(Tran21)./5);
    A21min = -(log(Tran21Min)./5);
    A21max = -(log(Tran21Max)./5);
    A43 = -(log(Tran43)./5);
    A43min = -(log(Tran43Min)./5);
    A43max = -(log(Tran43Max)./5);   
 
%     AEM = abs(log(TEM)./5);
%     AB1 = abs(log(TB1)./5);
%     AB0 = abs(log(TB0)./5);
%     A21 = abs(log(Tran21)./5);
%     A21min = abs(log(Tran21Min)./5);
%     A21max = abs(log(Tran21Max)./5);
%     A43 = abs(log(Tran43)./5);
%     A43min = abs(log(Tran43Min)./5);
%     A43max = abs(log(Tran43Max)./5); 
    
    
    %Figure - No Log Scale
    figure();
    p1 = plot(ffmodp,AEM,'--k', 'DisplayName','Model TDEMM');
    hold on;
    p2 = plot(ffmodp,AB1,'-.r', 'DisplayName','Model B0');
    p3 = plot(ffmodp,AB0,'-.b', 'DisplayName','Model B1');
    p4 = plot(ffp, A21, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
    ATWmaxp = A21max;
    ATWminp = A21min;
    ShadeX = [ffp(2:end); flipud(ffp(2:end))];
    inBetween = [(ATWmaxp(2:end)); flipud(ATWminp(2:end))];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    alpha(0.3)
    p4a = plot(ffp, ATWmaxp, ':r', 'LineWidth', 1, 'DisplayName','max/min');
    plot(ffp, ATWminp, ':r', 'LineWidth', 1);    
    
    p5 = plot(ffp, A43, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
    ATWmaxp = A43max;
    ATWminp = A43min;
    ShadeX = [ffp(2:end); flipud(ffp(2:end))];
    inBetween = [(ATWmaxp(2:end)); flipud(ATWminp(2:end))];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    alpha(0.3)
    p5a = plot(ffp, ATWmaxp, ':b', 'LineWidth', 1, 'DisplayName','max/min');
    plot(ffp, ATWminp, ':b', 'LineWidth', 1); 
    
    xlabel('Frequency / Fp');
    ylabel('Attenuation Coefficient ');
    title(['Attenuation  of Spectra Comparisons ', Description]);
    legend([p1 p2 p3 p4 p4a p5 p5a]);
    hold off; 
    
    %Log Sclae
    figure();
    p1 = plot(ffmodp,AEM,'--k', 'DisplayName','Model TDEMM');
    hold on;
    p2 = plot(ffmodp,AB1,'-.r', 'DisplayName','Model B0');
    p3 = plot(ffmodp,AB0,'-.b', 'DisplayName','Model B1');
    p4 = plot(ffp, A21, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
    ATWmaxp = A21max;
    ATWminp = A21min;
    ShadeX = [ffp(2:end); flipud(ffp(2:end))];
    inBetween = [(ATWmaxp(2:end)); flipud(ATWminp(2:end))];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    alpha(0.3)
    p4a = plot(ffp, ATWmaxp, ':r', 'LineWidth', 1, 'DisplayName','max/min');
    plot(ffp, ATWminp, ':r', 'LineWidth', 1);    
    p5 = plot(ffp, A43, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
    ATWmaxp = A43max;
    ATWminp = A43min;
    ShadeX = [ffp(2:end); flipud(ffp(2:end))];
    inBetween = [(ATWmaxp(2:end)); flipud(ATWminp(2:end))];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    alpha(0.3)
    p5a = plot(ffp, ATWmaxp, ':b', 'LineWidth', 1, 'DisplayName','max/min');
    plot(ffp, ATWminp, ':b', 'LineWidth', 1); 
    xlabel('Frequency / Fp');
    ylabel('Attenuation Coefficient ');
    title(['Attenuation  of Spectra Comparisons Log Scale', Description]);
    legend([p1 p2 p3 p4 p4a p5 p5a]);
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold off; 


    y12 = log(A21(and(and(ffp> 1, ffp < 2),A21 > 0)));
    x12 = log(ffp(and(and(ffp> 1, ffp < 2),A21 > 0)));    
    mdl12b12 = fitlm(x12,y12);  
    b12b12 = mdl12b12.Coefficients.Estimate(1);
    m12b12 = mdl12b12.Coefficients.Estimate(2);       
    
    y43 = log(A43(and(and(ffp> 1, ffp < 2),A43 > 0)));
    x43 = log(ffp(and(and(ffp> 1, ffp < 2),A43 > 0)));    
    mdl43b12 = fitlm(x43,y43);  
    b43b12 = mdl43b12.Coefficients.Estimate(1);
    m43b12 = mdl43b12.Coefficients.Estimate(2);    

    
    yEM = log(AEM(and(and(ffmodp> 0.8, ffmodp < 1.6),AEM > 0)));
    xEM = log(ffmodp(and(and(ffmodp> 0.8, ffmodp < 1.6),AEM > 0)));    
    mdlEMb12 = fitlm(xEM,yEM);  
    bEMb12 = mdlEMb12.Coefficients.Estimate(1);
    mEMb12 = mdlEMb12.Coefficients.Estimate(2);    
    
%     y12 = log(A21(and(and(ffp> 2, ffp < 3),A21 > 0)));
%     x12 = log(ffp(and(and(ffp> 2, ffp < 3),A21 > 0)));    
%     mdl12b23 = fitlm(x12,y12);  
%     b12b23 = mdl12b23.Coefficients.Estimate(1);
%     m12b23 = mdl12b23.Coefficients.Estimate(2);       
%     
%     y43 = log(A43(and(and(ffp> 2, ffp < 3),A43 > 0)));
%     x43 = log(ffp(and(and(ffp> 2, ffp < 3),A43 > 0)));    
%     mdl43b23 = fitlm(x43,y43);  
%     b43b23 = mdl43b23.Coefficients.Estimate(1);
%     m43b23 = mdl43b23.Coefficients.Estimate(2);    
    
    
    figure();
    p1 = plot(ffmodp,AEM,'--k', 'DisplayName','Model TDEMM');
    hold on;
    p2 = plot(ffmodp,AB1,'-.r', 'DisplayName','Model B0');
    p3 = plot(ffmodp,AB0,'-.b', 'DisplayName','Model B1');
    p4 = plot(ffp, abs(log(Tran21)/5), '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
    p5 = plot(ffp, abs(log(Tran43)/5), '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
%     p6 = plot(exp(x12), exp(y12), 'DisplayName','data interp 21');
%     p7 = plot(exp(x43), exp(y43), 'DisplayName','data interp 43');
    p8 = plot(ffp, exp(m12b12*log(ffp) + b12b12), 'DisplayName',['Model Slope (1,2) : ', num2str(m12b12)]);
    p9 = plot(ffp, exp(m43b12*log(ffp) + b43b12), 'DisplayName',['Model Slope (1,2): ', num2str(m43b12)]);
    p10 = plot(ffp, exp(mEMb12*log(ffp) + bEMb12), 'DisplayName',['Model Slope (EM) : ', num2str(mEMb12)]);
    xlabel('Frequency / Fp');
    ylabel('Attenuation Coefficient ');
    title(['Attenuation  of Spectra Comparisons Log Scale', Description]);
    legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10]);
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlim([0.5,2])
    hold off; 
 
    
    
    
    
    %Figure - Attenuation vs Incoming Amplitude - useless because our
    %amplitude is a function of frequency. 
%     figure();
%     plot(Amp1(Amp1 > Hs12/100.0 )./(4*sqrt(m01)), A21(Amp1 > Hs12/100.0 ), '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
%     hold on;
%     plot(Amp3(Amp3 > Hs34/100.0 )./(4*sqrt(m03)), A43(Amp3 > Hs34/100.0), '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
%     xlabel('Amplitude / Hs');
%     ylabel('Attenuation Coefficient ');
%     title(['Attenuation vs Amplitude ', Description]);
%     legend();
%     hold off; 

    
    %
%        
%     y12 = log(A21);
%     x12 = log(ffp);    
%     mdl12 = fitlm(x12,y12);  
%     b12 = mdl12.Coefficients.Estimate(1);
%     m12 = mdl12.Coefficients.Estimate(2);   
%     
%     y34 = log(A43);
%     x34 = log(ffp);    
%     mdl34 = fitlm(x34,y34);  
%     b34 = mdl34.Coefficients.Estimate(1);
%     m34 = mdl34.Coefficients.Estimate(2);     
% 
%     
%     figure();
%     plot(1./ffmodp,AEM,'--k', 'DisplayName','Model TDEMM');
%     hold on;
%     plot(1./ffmodp,AB1,'-.b', 'DisplayName','Model B0');
%     plot(1./ffmodp,AB0,'-.r', 'DisplayName','Model B1');
%     plot(1./ffp, A21, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
%     plot(1./ffp, A43, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
% %     plot(1./exp([-10;x12;10]), 1./exp(m12.*[-10;x12;10] + b12) ,'--r','DisplayName', ['LR12  Slope: ', num2str(m12)]);
% %     plot(1./exp([-10;x34;10]), 1./exp(m34.*[-10;x34;10] + b34) ,'--b','DisplayName', ['LR34  Slope: ', num2str(m34)]);
%     xlabel('T / Tp');
%     ylabel('Attenuation Coefficient ');
%     title(['Attenuation  of Spectra Comparisons ', Description]);
%     legend();
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
%     hold off;    

%     
%     y12 = log(A21ar1);
%     x12 = log(pff21ar1);    
%     mdl12 = fitlm(x12,y12);  
%     b12 = mdl12.Coefficients.Estimate(1);
%     m12 = mdl12.Coefficients.Estimate(2);   
% 
%     y34 = log(A43ar1);
%     x34 = log(pff43ar1);    
%     mdl34 = fitlm(x34,y34);  
%     b34 = mdl34.Coefficients.Estimate(1);
%     m34 = mdl34.Coefficients.Estimate(2);   
%     
%     figure();
%     plot(ffAEM,AEM,'--k', 'DisplayName','Model TDEMM');
%     hold on;
%     plot(ffAB1,AB1,'-.b', 'DisplayName','Model B0');
%     plot(ffAB0,AB0,'-.r', 'DisplayName','Model B1');
%     plot(pff21, A21, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
%     plot(pff21p,A21p,'-r')
%     plot(pff21m,A21m,'-r')
%     plot(pff43, A43, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
%     plot(pff43p,A43p,'-b')
%     plot(pff43m,A43m,'-b')
%     plot(exp([-10;x12;10]), exp(m12.*[-10;x12;10] + b12) ,'--r','DisplayName', ['LR12  Slope: ', num2str(m12)]);
%     plot(exp([-10;x34;10]), exp(m34.*[-10;x34;10] + b34) ,'--b','DisplayName', ['LR34  Slope: ', num2str(m34)]);
%     %plot(exp([-10;x12;10]), exp(3.*[-10;x12;10] + b12) ,'--k','DisplayName', ['Target  Slope 3 Over 12 ']);
%     %plot(exp([-10;x34;10]), exp(3.*[-10;x34;10] + b34) ,'--k','DisplayName', ['Target  Slope 3 Over 34 ']);
%     plot(pff21ar1,A21ar1);
%     plot(pff43ar1,A43ar1);
%     xlabel('Frequency / Fp');
%     ylabel('Attenuation Coefficient ');
%     title(['Attenuation  of Spectra Comparisons ', Description]);
%     legend();
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
%     hold off;    

    
%     %data
%     meanalpha12 = mean(alpha12(and(ff1{1} >Tp/2,ff1{1} < 6*Tp/4 )));
%     meanalpha34 = mean(alpha34(and(ff1{1} >Tp/2,ff1{1} < 6*Tp/4 )));
%     MeanValueAtTp = mean([meanalpha12,meanalpha34]);
%     
%     %linear regression
%     y12 = log(alpha12(and(ff12 >Tp/2,ff12 < 6*Tp/4 )));
%     x12 = log(ff12(and(ff12  >Tp/2,ff12 < 6*Tp/4 ))/Tp);
%     mdl12 = fitlm(y12,x12);
%     m12 = mdl12.Coefficients.Estimate(1);
%     b12 = mdl12.Coefficients.Estimate(2);
%     
%     %linear regression
%     y34 = log(alpha34(and(ff34 >Tp/2,ff34 < 6*Tp/4 )));
%     x34 = log(ff34(and(ff34  >Tp/2,ff34 < 6*Tp/4 ))/Tp);
%     mdl34 = fitlm(y34,x34);
%     m34 = mdl34.Coefficients.Estimate(1);
%     b34 = mdl34.Coefficients.Estimate(2);
%     
%     %Model
%     
%     %linear regression
%     yTDEMM = log(alphaTDEMM);
%     yB1 = log(alphaB1);
%     yB0 = log(alphaB0);
%     xff = log(ffmodl);
%     
%     mdlTDEMM = fitlm(yTDEMM,xff);
%     mdlB0 = fitlm(yB0,xff);
%     mdlB1 = fitlm(yB1,xff);
% 
%     figure();
%     plot(ff12/Tp, alpha12, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
%     hold on;
%     plot(ff34/Tp, alpha34, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
% %     plot(ff1{1}/Tp, (MeanValueAtTp)*(ff1{1}/Tp).^3 , '--k', 'DisplayName',  ['Power Law Alpha = 3'] , 'LineWidth', 2);
% %     plot(ff1{1}/Tp, (MeanValueAtTp)*(ff1{1}/Tp).^2 , '--g', 'DisplayName',  ['Power Law Alpha = 2'] , 'LineWidth', 2);
% %     plot(ff1{1}/Tp, (MeanValueAtTp)*(ff1{1}/Tp).^4 , ':k', 'DisplayName',  ['Power Law Alpha = 4'] , 'LineWidth', 2);
% 
%     %     plot(ff1{1}*Tp, (ff1{1}*Tp).^3, '--k', 'DisplayName',  ['Power Law Alpha = 3'] , 'LineWidth', 2);
% %     plot(ff1{1}*Tp, (ff1{1}*Tp).^2, '--g', 'DisplayName',  ['Power Law Alpha = 2'] , 'LineWidth', 2);
%     xlabel('Frequency');
%     ylabel('Attenuation ');
%     title(['Attenuation of Spectra Comparisons ', Description]);
%     legend();
% %     set(gca, 'YScale', 'log')
% %     set(gca, 'XScale', 'log')
%     xlim([0,2])
%     hold off;
% %     ylim([10^-12 , 1 ]);
% 
%     figure();
%     plot(xff,yTDEMM,'-r');
%     hold on;
%     plot(xff,yB0,'-b');
%     plot(xff,yB1,'-k');
%     plot(xff,xff,'--k');
%     plot(ff12*Tp, alpha12, '.r', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(Hs12),' [m]'] , 'MarkerSize', 12);
%     plot(ff34*Tp, alpha34, '.b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Hs ', num2str(Hs34),' [m]'] , 'MarkerSize', 12);
%     
%     
%     
%     figure();
%     plot(xff,yTDEMM,'-r');
%     hold on;
%     plot(xff,yB1,'-b');
%     plot(xff,yB0,'-k');
% %     plot(ff1{1}*Tp, ATW3./m03, '-b', 'DisplayName',  ['Conc ', num2str(conc{3}) ,' | Probes : ( ', num2str(probes{3}(1)),' , ', num2str(probes{3}(end)) ' ) | Hs ', num2str(Hs34),' [m]'] , 'LineWidth', 2);
% %     plot(ffJS*Tp, JSspec12./(Hs12/4)^2, '--k', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs12),' [m]'] , 'LineWidth', 2);
% %     plot(ffJS*Tp, JSspec34./(Hs34/4)^2, '--g', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs34),' [m]'] , 'LineWidth', 2);
% %     xlabel('Frequency / Fp');
% %     ylabel('Spectra / m0 ');
% %     title(['Calibration Spectra Comparisons', Description]);
% %     legend();
% %     set(gca, 'YScale', 'log')
% %     ylim([10^-12 , 1 ]);
% %     hold off; 
%      
%     
%     figure();
%     for i = 1: length( SnS1)
%         if i == 1
%             p1l = plot(ff1{i}*Tp, SnS1{i}, '-r', 'DisplayName', ['Conc ', num2str(conc{1})  ' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' )']);
%             hold on;
%             p2l = plot(ff2{i}*Tp, SnS2{i}, '-b', 'DisplayName', ['Conc ', num2str(conc{2})  ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )']);
%         end
%             plot(ff1{i}*Tp, SnS1{i}, '-r')
%             plot(ff2{i}*Tp, SnS2{i}, '-b')
%     end
%     xlabel('Frequency / Fp');
%     ylabel('Spectra / Hs');
%     title(['Individual Probes (Average over Time Intervals) ', Description]);
%     legend([p1l,p2l]);
%     set(gca, 'YScale', 'log')
%     ylim([10^-12 , 10^-2 ]);
%     hold off
%     
%     %Plot of spectra and predictions
%     
%     %Compare Calib and Conc


%     figure();
%     plot(ffJS, JSspec, '-r', 'DisplayName',  'JonSwap' , 'LineWidth', 2);
%     hold on;
%     plot(ff1{1}, AverageTimeWindows2, '-b', 'DisplayName',  ['Conc ', num2str(conc{2}), ' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' )'] , 'LineWidth', 2)    ;
%     plot(ffJS, PTC_2EEMJS.*JSspec, '--k', 'DisplayName', '2DEMM Predicted Spectra' , 'LineWidth', 2);
%     plot(ffJS, PTC_Bolt0JS.*JSspec, '-.r', 'DisplayName', 'Boltzman Predicted Spectra' , 'LineWidth', 2);
%     plot(ffJS, PTC_Bolt1JS.*JSspec, '-.g', 'DisplayName', 'Boltzman (Remove Scattering) Predicted Spectra' , 'LineWidth', 2);
%     xlabel('Frequency [Hz]');
%     ylabel('Spectra ');
%     title(['Average Spectra (Probes and Time Intervals) with Predicted JONSWAP Spectra ', Description]);
%     legend();
%     xlim([0 1.0/ (Tp/4)]);
%     hold off;
%     
%     %Transmission Coefficients
%     figure()
%     plot(1./ff1{1}, AverageTimeWindows2 ./ AverageTimeWindows1,'xr','DisplayName','Data' );
%     hold on;
%     plot(model_pers,TwoDEMM,'--k','DisplayName','2D EEM' )
%     plot(model_pers,Boltzman0,'--b','DisplayName','Boltzman' )
%     plot(model_pers,Boltzman1,'--r','DisplayName','Boltzman (Remove Scattering)' )
%     xlabel('Periods');
%     ylabel('Transmission Coefficient');
%     title(['Transmission Coefficient', Description]);
%     legend();
%     xlim([0,3]);  
%     ylim([0,1]);
%     hold off;
%     
%     figure()
%     plot(1./ff1{1}, AverageTimeWindows2 ./ JSspec1,'xr','DisplayName','Data' );
%     hold on;
%     plot(model_pers,TwoDEMM,'--k','DisplayName','2D EEM' )
%     plot(model_pers,Boltzman0,'--b','DisplayName','Boltzman' )
%     plot(model_pers,Boltzman1,'--r','DisplayName','Boltzman (Remove Scattering)' )
%     xlabel('Periods (s)');
%     ylabel('Transmission Coefficient');
%     title(['Transmission Coefficient JonSwap', Description]);
%     legend();
%     xlim([0,3]);  
%     ylim([0,1]);
%     hold off;
%     
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

