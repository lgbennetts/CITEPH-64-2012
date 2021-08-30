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

function fn_ManyProbes(conc,TestName,Tp,probes)

%Compare calibration to disk
if ~exist('conc','var') conc={0,0}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'4','4'}}; end ;  % 9, 14 %14 %2

% if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  

% if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  % 9, 14


if ~exist('Cols','var'); Cols= {'xr', '^b','sk'}; end ;

if ~exist('wbinwidth','var');  wbinwidth=6; end;

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
        
        Hs = c_pram.wave_height/100.0;

        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType);
        t0index = find(strcmp({Tind.description},t0description));
        t0 = Tind(t0index).time + T_windowSec;
        t1index = find(strcmp({Tind.description},t1description));
        t1 = Tind(t1index).time - T_windowSec;

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp); % fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        jj0=find(tm_vec>t0,1);
        jj1=find(tm_vec<t1,1,'last');
        
        t0jj0 = tm_vec(jj0);
        t1jj1 = tm_vec(jj1);
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        SnS1{i} = Sn1; 
        
        %Raw signal figure
%         if (i == 1)
%             figure();
%             plot(tm,disp,'-k')
%             hold on;
%             p1 = plot([t0jj0 t0jj0], [min(disp), max(disp)], '--b' , 'DisplayName', ['Time Interval Begin']);
%             p2 = plot([t1jj1 t1jj1], [min(disp), max(disp)], '--r' , 'DisplayName', ['Time Interval End']);
%             legend([p1 p2])
%             ylabel('Displacement [m]')
%             xlabel('Time [s]')
%             title(['Raw Signal Hs: ', num2str(Hs), ' Tp :', num2str(Tp), ' Probe ', num2str(probes{1}(i))  ])
%             hold off;
%         end

        %% Fourier Transform


        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{j}{2}, probes{2}(i));

        Tp = c_pram.period/10.0;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType);
        t0index = find(strcmp({Tind.description},t0description));
        t0 = Tind(t0index).time + T_windowSec;
        t1index = find(strcmp({Tind.description},t1description));
        t1 = Tind(t1index).time - T_windowSec;
        
        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp); %fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;
        
        jj0=find(tm_vec>t0,1);
        jj1=find(tm_vec<t1,1,'last');
        
        Sn2 = mean(Sn_mat(:,jj0:jj1), 2);
        SnS2{i} = Sn2;

        i = i +1;
    end

%     ffp = ff*Tp;
    ffpA = ff.*Tp;
    ffA = ff;
    SnS1L = [];
    SnS2L = [];
    for i = 1: length( SnS1)
        SnS1L = [SnS1L , SnS1{i} ];
        SnS2L = [SnS2L , SnS2{i} ];
    end
   
    ATWR = mean(SnS1L,2);
    ATWT = mean(SnS2L,2);
    
    ATWRmax = max(SnS1L,[],2);
    ATWRmin = min(SnS1L,[],2);
    ATWTmax = max(SnS2L,[],2);
    ATWTmin = min(SnS2L,[],2);

    wJS = 2*pi/Tp;
    ffJS = ff;
    JSspec  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
    ffJSp = ffJS*Tp;
    
    Description = ['Concentations (', num2str(conc{1}) ,' , ', num2str(conc{2}) ') Tp = ', num2str(Tp), '[s] Hs = ', num2str(Hs) , '[m] Fourier Interval ', num2str(FourierWindow_Sec),'[s] ' ];

    
    m0R = trapz(ffA,ATWR);
    m0T = trapz(ffA,ATWT);
    
    
    %Smooth spectra
    wJS = 2*pi/Tp;
    ffJSsmoothnn = 0:0.01:100;
    JSSmooth  = jonswap(2*pi*ffJSsmoothnn,'wp',wJS,'Hs',Hs);
    ffJSSmooth = ffJSsmoothnn*Tp;
    
    PlotNames = {};
    %Plot Mean Spectra with Maximum and Minimum
    figure();
    set(gca,'FontSize',18) 
    p1 = plot(ffpA, ATWR./m0R, '-r' , 'LineWidth', 2);   
    PlotNames{end +1} = ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(Hs),' [m]'];
    hold on;
    ATWRmaxplot = ATWRmax(2:end)./m0R;
    ATWRminplot = ATWRmin(2:end)./m0R;
    ShadeX = [ffpA(2:end); flipud(ffpA(2:end))];
    inBetween = [ATWRmaxplot; flipud(ATWRminplot)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none')
    PlotNames{end +1} = 'Max and Min Region (Among Probes)';
    alpha(0.3) 
    p1a = plot(ffpA(2:end), ATWRmaxplot, ':r' , 'LineWidth', 2);
    PlotNames{end +1} = 'Maximum and Minimum';
    
    p1b = plot(ffpA(2:end), ATWRminplot, ':r' , 'LineWidth', 2);
    PlotNames{end +1} = '';
    
    p2 = plot(ffpA, ATWT./m0T, '-b' , 'LineWidth', 2);
    PlotNames{end +1} = ['Conc ', num2str(conc{2}) ,' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(Hs),' [m]'];
    ATWTmaxplot = ATWTmax(2:end)./m0T;
    ATWTminplot = ATWTmin(2:end)./m0T;
    ShadeX = [ffpA(2:end); flipud(ffpA(2:end))];
    inBetween = [ATWTmaxplot; flipud(ATWTminplot)];
    patch(ShadeX, inBetween, 'b', 'LineStyle','none')
    PlotNames{end +1} = 'Max and Min Region (Among Probes)';
    alpha(0.3) 
    p2a = plot(ffpA(2:end), ATWTmaxplot, ':b', 'DisplayName', 'Maximum and Minimum' , 'LineWidth', 2);
    PlotNames{end +1} = 'Maximum and Minimum';
    p2b = plot(ffpA(2:end), ATWTminplot, ':b' , 'LineWidth', 2);
    PlotNames{end +1} = '';
    
    xlabel('Frequency / Fp');
    ylabel('S(f) / m0 ');
    xlim([0,3])
    title(['Spectra Comparison of Probes', Description]);
    legend(PlotNames);
    hold off; 
    
    HserrorR = abs(4*sqrt(m0R) - Hs) / Hs
    HserrorT = abs(4*sqrt(m0T) - Hs) / Hs

    a1 = gca;
    f2 = figure;
    a2 = copyobj(a1,f2) ;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlim([10^-1 , 10^2 ]);
    ylim([10^-8 , 10 ]);
    legend(PlotNames);
end

%     
%     Mean Spectra with min and max
%     figure();
%     p1 = plot(ffpA(:,1), ATWR(:,1)./m0R(1,1), '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA(1,1)),' [m]'] , 'LineWidth', 3);
%     hold on;
%    
% 
%     
%     p2 = plot(ffpA(:,2), ATWR(:,2)./m0R(1,2), '-b', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA(1,2)),' [m]'] , 'LineWidth', 2);
%     
%     ATWRmaxplot = ATWRmax(:,2)./m0R(1,2);
%     ATWRminplot = ATWRmin(:,2)./m0R(1,2);
%     ShadeX = [ffpA(:,2); flipud(ffpA(:,2))];
%     inBetween = [ATWRmaxplot; flipud(ATWRminplot)];
%     patch(ShadeX, inBetween, 'b', 'LineStyle','none')
%     alpha(0.3) 
%     p2a = plot(ffpA(:,2), ATWRmaxplot, ':b', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p2b = plot(ffpA(:,2), ATWRminplot, ':b', 'DisplayName', 'Min' , 'LineWidth', 2);
%     
%     p3 = plot(ffpA(:,1), ATWT(:,1)./m0T(1,1), '-','color','#ff8d74', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(HsA(1,1)),' [m]'] , 'LineWidth', 2);
%     
%     ATWTmaxplot = ATWTmax(:,1)./m0T(1,1);
%     ATWTminplot = ATWTmin(:,1)./m0T(1,1);
%     ShadeX = [ffpA(:,1); flipud(ffpA(:,1))];
%     inBetween = [ATWTmaxplot; flipud(ATWTminplot)];
%     patch(ShadeX, inBetween,'r','FaceColor','#ff8d74', 'LineStyle','none')
%     alpha(0.3)  
%     p3a = plot(ffpA(:,1), ATWTmaxplot, ':','color','#ff8d74', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p3b = plot(ffpA(:,1), ATWTminplot, ':','color','#ff8d74', 'DisplayName', 'Min' , 'LineWidth', 2);
%     
%     p4 = plot(ffpA(:,2), ATWT(:,2)./m0T(1,2), '-','color','#74e2ff', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(HsA(1,2)),' [m]'] , 'LineWidth', 2);
% 
%     ATWTmaxplot = ATWTmax(:,2)./m0T(1,2);
%     ATWTminplot = ATWTmin(:,2)./m0T(1,2);
%     ShadeX = [ffpA(:,2); flipud(ffpA(:,2))];
%     inBetween = [ATWTmaxplot; flipud(ATWTminplot)];
%     patch(ShadeX, inBetween,'b','FaceColor','#74e2ff', 'LineStyle','none')
%     alpha(0.3)  
%     p4a = plot(ffpA(:,2), ATWTmaxplot, ':','color','#74e2ff', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p4b = plot(ffpA(:,2), ATWTminplot, ':','color','#74e2ff', 'DisplayName', 'Min' , 'LineWidth', 2);
%     
%     p5 = plot(ffJSp(:,j), JSspec(:,j)./(Hs/4)^2, '--k', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs),' [m]'] , 'LineWidth', 2);
%     xlabel('Frequency / Fp');
%     ylabel('Spectra / m0 ');
%     title(['Calibration Spectra Comparisons With Max and Min', Description]);
%     legend([p1 p1a p1b p2 p2a p2b p3 p3a p3b p4 p4a p4b p5]);
%     xlim([0.5 , 2]);
%     hold off;    
% 
%     Mean Spectra with max and min
%     figure();
%     p1 = plot(ffpA(:,1), ATWR(:,1)./m0R(1,1), '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA(1,1)),' [m]'] , 'LineWidth', 3);
%     hold on;
%     
%     ATWRmaxplot = ATWRmax(2:end,1)./m0R(1,1);
%     ATWRminplot = ATWRmin(2:end,1)./m0R(1,1);
%     ShadeX = [ffpA(2:end,1); flipud(ffpA(2:end,1))];
%     inBetween = [ATWRmaxplot; flipud(ATWRminplot)];
%     patch(ShadeX, inBetween, 'r', 'LineStyle','none')
%     alpha(0.3)
%     p1a = plot(ffpA(2:end,1), ATWRmaxplot, ':r', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p1b = plot(ffpA(2:end,1), ATWRminplot, ':r', 'DisplayName', 'Min' , 'LineWidth', 2);
% 
%     
%     p2 = plot(ffpA(:,2), ATWR(:,2)./m0R(1,2), '-b', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA(1,2)),' [m]'] , 'LineWidth', 2);
%     
%     ATWRmaxplot = ATWRmax(2:end,2)./m0R(1,2);
%     ATWRminplot = ATWRmin(2:end,2)./m0R(1,2);
%     ShadeX = [ffpA(2:end,2); flipud(ffpA(2:end,2))];
%     inBetween = [ATWRmaxplot; flipud(ATWRminplot)];
%     patch(ShadeX, inBetween, 'b', 'LineStyle','none')
%     alpha(0.3) 
%     p2a = plot(ffpA(2:end,2), ATWRmaxplot, ':b', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p2b = plot(ffpA(2:end,2), ATWRminplot, ':b', 'DisplayName', 'Min' , 'LineWidth', 2);
%     
%     p3 = plot(ffpA(:,1), ATWT(:,1)./m0T(1,1), '-','color','#ff8d74', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(HsA(1,1)),' [m]'] , 'LineWidth', 2);
%     
%     ATWTmaxplot = ATWTmax(2:end,1)./m0T(1,1);
%     ATWTminplot = ATWTmin(2:end,1)./m0T(1,1);
%     ShadeX = [ffpA(2:end,1); flipud(ffpA(2:end,1))];
%     inBetween = [ATWTmaxplot; flipud(ATWTminplot)];
%     patch(ShadeX, inBetween,'r','FaceColor','#ff8d74', 'LineStyle','none')
%     alpha(0.3)  
%     p3a = plot(ffpA(2:end,1), ATWTmaxplot, ':','color','#ff8d74', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p3b = plot(ffpA(2:end,1), ATWTminplot, ':','color','#ff8d74', 'DisplayName', 'Min' , 'LineWidth', 2);
%     
%     p4 = plot(ffpA(:,2), ATWT(:,2)./m0T(1,2), '-','color','#74e2ff', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(HsA(1,2)),' [m]'] , 'LineWidth', 2);
% 
%     ATWTmaxplot = ATWTmax(2:end,2)./m0T(1,2);
%     ATWTminplot = ATWTmin(2:end,2)./m0T(1,2);
%     ShadeX = [ffpA(2:end,2); flipud(ffpA(2:end,2))];
%     inBetween = [ATWTmaxplot; flipud(ATWTminplot)];
%     patch(ShadeX, inBetween,'b','FaceColor','#74e2ff', 'LineStyle','none')
%     alpha(0.3)  
%     p4a = plot(ffpA(2:end,2), ATWTmaxplot, ':','color','#74e2ff', 'DisplayName', 'Max' , 'LineWidth', 2);
%     p4b = plot(ffpA(2:end,2), ATWTminplot, ':','color','#74e2ff', 'DisplayName', 'Min' , 'LineWidth', 2);
%     
%     p5 = plot(ffJSp(:,j), JSspec(:,j)./(Hs/4)^2, '--k', 'DisplayName',  ['JonSwap | Hs ', num2str(Hs),' [m]'] , 'LineWidth', 2);
%     xlabel('Frequency / Fp');
%     ylabel('Spectra / m0 ');
%     title(['Calibration Spectra Comparisons Log Scales', Description]);
%     legend([p1 p1a p1b p2 p2a p2b p3 p3a p3b p4 p4a p4b p5]);
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
%     ylim([10^-12 , 10]);
%     hold off;   
%     
% 
%     PLot of transmission / reflection
%     figure();
%     p1 = plot(ffpA(:,1), abs(ATWT(:,1) - ATWR(:,1)) ./ ATWR(:,1), '-r', 'DisplayName',  ['Conc ', num2str(conc{1}) ,' | Hs ', num2str(HsA(1,1)),' [m]'] , 'LineWidth', 3);
%     hold on;
%     p2 = plot(ffpA(:,2), abs(ATWT(:,2)- ATWR(:,2)) ./ ATWR(:,2), '-b', 'DisplayName',  ['Conc ', num2str(conc{2}) ,' | Hs ', num2str(HsA(1,2)),' [m]'] , 'LineWidth', 3);
%     xlabel('Frequency / Fp');
%     ylabel('Transmission / Reflection ');
%     title(['Relative Error Spectra (Reflection and Transmission) Comparisons Log Scales ', Description]);
%     legend([p1 p2]);
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
%     ylim([10^-4 , 100]);
%     hold off;   
% 


return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

