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

function fn_Transmission(conc,TestName,probes)

% %Compare calibration to disk
% if ~exist('conc','var') conc={0,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'16','8'}}; end ;  
% if ~exist('TestName','var'); TestName= {{'14','9'},{'15','10'},{'16','8'}}; end ;  
% 
if ~exist('conc','var') conc={0,79}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'16','21'}}; end ;  
% if ~exist('TestName','var'); TestName= {{'14','19'}, {'15','20'},{'16','21'},{'17','22'}}; end ;  

% if ~exist('Cols','var'); Cols= {'xr', '^b','sk','dg'}; end ;
if ~exist('Cols','var'); Cols= {{'#ff0000','#680303'},{'#0bff01 ','#058000'},{'#0487f9','#00427c'},{'#9701ff','#470179'}}; end ;


if ~exist('binnum','var');  binnum=30; end;

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
        %Read Data in Calibration Reflection Probes
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, probes{1}(i));

        Tp = c_pram.period/10.0;
        
        dt = tm(2) - tm(1);   
        
        eval(TpersStr);
        T_windowSec = round(Tpers*Tp);
        SamplingRate = floor(1./dt);
        T_window = SamplingRate*T_windowSec;
        
        Hs = c_pram.wave_height/100.0;
        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        t0index = find(strcmp({Tind.description},t0description));
        t0 = Tind(t0index).time + T_windowSec;
        t1index = find(strcmp({Tind.description},t1description));
        t1 = Tind(t1index).time - T_windowSec ;

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        jj0=find(tm_vec>t0,1);
        jj1=find(tm_vec<t1,1,'last');
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        
        ffCalibPR{i} = ff;
        SnSCalibPR{i} = Sn1;
        TimesCalibPR{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];        
        

        %Read Data in Calibration Transmission Probes
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, probes{2}(i));

        Tp = c_pram.period/10.0;
        
        dt = tm(2) - tm(1);   
        
        eval(TpersStr);
        T_windowSec = round(Tpers*Tp);
        SamplingRate = floor(1./dt);
        T_window = SamplingRate*T_windowSec;
        
        Hs = c_pram.wave_height/100.0;
        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        t0index = find(strcmp({Tind.description},t0description));
        t0 = Tind(t0index).time + T_windowSec;
        t1index = find(strcmp({Tind.description},t1description));
        t1 = Tind(t1index).time - T_windowSec ;

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        jj0=find(tm_vec>t0,1);
        jj1=find(tm_vec<t1,1,'last');
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        
        ffCalibPT{i} = ff;
        SnSCalibPT{i} = Sn1;
        TimesCalibPT{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];  
        

        
         %Read Data in Experiment Reflection Probes
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{j}{2}, probes{1}(i));

        Tp = c_pram.period/10.0;
        
        dt = tm(2) - tm(1);   
        
        eval(TpersStr);
        T_windowSec = round(Tpers*Tp);
        SamplingRate = floor(1./dt);
        T_window = SamplingRate*T_windowSec;
        
        Hs = c_pram.wave_height/100.0;
        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        t0index = find(strcmp({Tind.description},t0description));
        t0 = Tind(t0index).time + T_windowSec;
        t1index = find(strcmp({Tind.description},t1description));
        t1 = Tind(t1index).time - T_windowSec ;

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        jj0=find(tm_vec>t0,1);
        jj1=find(tm_vec<t1,1,'last');
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        
        ffExpPR{i} = ff;
        SnSExpPR{i} = Sn1;
        TimesExpPR{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];        
        

        %Read Data in Calibration Transmission Probes
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{j}{2}, probes{2}(i));

        Tp = c_pram.period/10.0;
        
        dt = tm(2) - tm(1);   
        
        eval(TpersStr);
        T_windowSec = round(Tpers*Tp);
        SamplingRate = floor(1./dt);
        T_window = SamplingRate*T_windowSec;
        
        Hs = c_pram.wave_height/100.0;
        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        t0index = find(strcmp({Tind.description},t0description));
        t0 = Tind(t0index).time + T_windowSec;
        t1index = find(strcmp({Tind.description},t1description));
        t1 = Tind(t1index).time - T_windowSec ;

        [Sn_mat,tm_vec,ff,FourierWindow] = fn_JustFourierMovingWindow(tm,disp,T_window,SamplingRate,Tp);
        period = 1./ff;
        FourierWindow_Sec = FourierWindow/SamplingRate;

        jj0=find(tm_vec>t0,1);
        jj1=find(tm_vec<t1,1,'last');
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        
        ffExpPT{i} = ff;
        SnSExpPT{i} = Sn1;
        TimesExpPT{i} = [round(tm_vec(jj0));round(tm_vec(jj1))];  

        i = i +1;
    end

    SnSCalibPTL = [];
    SnSCalibPRL = [];
    SnSExpPTL = [];
    SnSExpPRL = [];
    for i = 1: length( SnSCalibPT)
        SnSCalibPTL = [SnSCalibPTL , SnSCalibPT{i} ];
        SnSCalibPRL = [SnSCalibPRL , SnSCalibPR{i} ];
        SnSExpPTL = [SnSExpPTL , SnSExpPT{i} ];
        SnSExpPRL = [SnSExpPRL , SnSExpPR{i} ];
    end

    %experimental
    TpA{j} = Tp;
    HsA{j} = Hs;
    ffA{j} = ff;
    ffpA{j} =ff*Tp;
   
    ATWCalibPR{j} = mean(SnSCalibPRL,2);
    ATWCalibPT{j} = mean(SnSCalibPTL,2);
    ATWExpPR{j} = mean(SnSExpPRL,2);
    ATWExpPT{j} = mean(SnSExpPTL,2);
    
%     
%     ATWCalibmax{j} = max(SnS1L,[],2);
%     ATWCalibmin{j} = min(SnS1L,[],2);
%     ATWConcmax{j} = max(SnS2L,[],2);
%     ATWConcmin{j} = min(SnS2L,[],2);
    
    RefEnergy{j} = max((ATWExpPR{j} - ATWCalibPR{j}) ./ ATWCalibPR{j} ,0);
    
    TransEnergy{j} = (ATWExpPT{j}) ./ ATWCalibPT{j}  ;
    
    InPlaneEnergy{j} =  RefEnergy{j} + TransEnergy{j};
    
%     fpL = 0;
%     fpH = 2;
%     AmmTol = 10^-3;
%     STol = AmmTol^2/2;
%     STol = 0;
    
    %above 1mm accuracy
%     Cond1 = ATWConc{j} > STol;
%     Cond2 = and(ffpA{j} > fpL, ffpA{j} < fpH);
%    
%     CondInd = and(Cond1,Cond2);
%     
%     ffFilt{j} = ffA{j}(CondInd);
%     ffpFilt{j} = ffpA{j}(CondInd);
%     ATWCalibFilt{j} = ATWCalib{j}(CondInd);    
%     ATWConcFilt{j} = ATWConc{j}(CondInd);
%     
%     ATWCalibFiltmax{j} = ATWCalibmax{j}(CondInd);    
%     ATWConcFiltmax{j} = ATWConcmax{j}(CondInd);
%     
%     ATWCalibFiltmin{j} = ATWCalibmin{j}(CondInd);    
%     ATWConcFiltmin{j} = ATWConcmin{j}(CondInd);
%     
%     Ref{j} = (ATWConcFilt{j}./ATWCalibFilt{j}) - 1;
%     RefMin{j} = (ATWConcFiltmin{j}./ATWCalibFilt{j}) - 1;
%     RefMax{j} = (ATWConcFiltmax{j}./ATWCalibFilt{j}) -1;
    
%     Alpha{j} = -log(Tran{j})./5;
    
    %target spectra
    wJS = 2*pi/Tp;
    ffJS = ff;
    JSspec{j}  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
    ffJSp{j} = ffJS*Tp;
    
    j = j +1;
end


Description = ['Concentration ' num2str(conc{2}) '% '];

%frequency plot range - smooth
ffJSsmoothnn = 0:0.01:100;

%Plot of in plane energy
figure();
set(gca,'FontSize',18) 
hold on;
PlotNames = {};
for j = 1: length(TestName)
    
    m0 = trapz(ffA{j},ATWCalibPR{j});
    
    plot(ffA{j}, InPlaneEnergy{j} , '-', 'Color',Cols{j}{1} , 'LineWidth', 3);
    pnameCalib = ['In Plane Energy | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
    PlotNames{end +1} = pnameCalib;
    
    plot(ffA{j}, ATWCalibPR{j}/m0, '-.b' , 'LineWidth', 3);
    pnameCalib = ['In Plane Energy | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
    PlotNames{end +1} = pnameCalib;    
end
plot([0,2], [1,1], '--k' , 'LineWidth', 3);
PlotNames{end +1} = 'Full In Plane';

xlabel('Frequency');
ylabel('Measured In Plane Energy');
title(['In Plane Energy', Description]);
legend(PlotNames);
                          
%Plot of Calibration vs Experiment Comparisons

%Spectra Comparisons - both Calibration and Experimental

for j = 1: length(TestName)
    figure();
    set(gca,'FontSize',18) 
    hold on;
    PlotNames = {};
    
    m0Calib = trapz(ffA{j},ATWCalib{j});
    
    plot(ffpA{j}, ATWCalib{j}./m0Calib , '-', 'Color',Cols{j}{1} , 'LineWidth', 3);
    pnameCalib = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
    PlotNames{end +1} = pnameCalib;
    
    Upper = ATWCalibmax{j}(2:end)./m0Calib;
    Lower = ATWCalibmin{j}(2:end)./m0Calib;
    ffrem0 = ffpA{j}(2:end);
    ShadeX = [ffrem0; flipud(ffrem0)];
    inBetween = [Upper ; flipud(Lower)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
    alpha(0.5)
    PlotNames{end +1} = 'Max and Min Region (Among Probes)';
    
    plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
    PlotNames{end +1} = 'Max and Min Boundary';
    plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
    PlotNames{end +1} = '';

    plot(ffpA{j}, ATWConc{j}./m0Calib , '-', 'Color',Cols{j}{2} , 'LineWidth', 3);
    pnameConc = ['Conc ', num2str(conc{2}),' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
    PlotNames{end +1} = pnameConc;
    
    Upper = ATWConcmax{j}(2:end)./m0Calib;
    Lower = ATWConcmin{j}(2:end)./m0Calib;
    ffrem0 = ffpA{j}(2:end);
    ShadeX = [ffrem0; flipud(ffrem0)];
    inBetween = [Upper ; flipud(Lower)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
    alpha(0.5)
    PlotNames{end +1} = 'Max and Min Region (Among Probes)';
    
    plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
    PlotNames{end +1} = 'Max and Min Boundary';
    plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
    PlotNames{end +1} = '';
    
    wJS = 2*pi/TpA{j};
    JSSmooth  = jonswap(2*pi*ffJSsmoothnn,'wp',wJS,'Hs',HsA{j});
    JSSmoothNorm = JSSmooth / (HsA{j}/4)^2;
    ffJSSmooth = ffJSsmoothnn*TpA{j};
    
    plot(ffJSSmooth, JSSmoothNorm, '--k' , 'LineWidth', 3);
    PlotNames{end +1} = 'JonSwap [m]';

    xlabel('Frequency / Fp');
    ylabel('Spectra / m0 ');
    title(['Calibration Spectra Comparisons', Description]);
    legend(PlotNames);
    xlim([0 , 3 ]);

    a1 = gca;
    f2 = figure;
    a2 = copyobj(a1,f2) ;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlim([10^-1 , 10^2 ]);
    ylim([10^-8 , 10 ]);
    legend(PlotNames);
end




% Models
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

ffmod = 1.0./ model_pers.';
ffmod = ffmod(2:end-1);
TEM = TwoDEMM(2:end-1);
TB0 = Boltzman0(2:end-1);
TB1 = Boltzman1(2:end-1);
AEM = -(log(TEM)./5);
AB1 = -(log(TB1)./5);
AB0 = -(log(TB0)./5);

%Transmission PLots
%Spectra Comparisons - both Calibration and Experimental
figure();
hold on;
 set(gca,'FontSize',18) 
clear PlotNames;
PlotNames = {};
for j = 1: length(TestName)
    plot(ffFilt{j}, Ref{j} , '.', 'Color',Cols{j}{1} ,'MarkerSize', 16);
    pname = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
    PlotNames{end +1} = pname;
    
    Upper = RefMax{j};
    Lower = RefMin{j};
    ffrem0 = ffFilt{j};
    ShadeX = [ffrem0; flipud(ffrem0)];
    inBetween = [Upper ; flipud(Lower)];
    patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
    alpha(0.3)
    PlotNames{end +1} = 'Max and Min Region (Among Probes)';
    
    plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
    PlotNames{end +1} = 'Max and Min Boundary';
    plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
    PlotNames{end +1} = '';
    
end
plot(ffmod, 1 -TB0 , '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(ffmod, 1-  TB1, '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(ffmod, 1 - TEM, '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

xlabel('Frequency');
ylabel('Reflection Coefficient ');
title(['Reflection Coefficient (Transmission Coefficient Model Overlaid) ', Description]);
legend(PlotNames);
xlim([0 , 5 ]);
ylim([0 , 1 ]);


%Attenuation Coefficient
%What if we just fit one?
figure();
hold on;
set(gca,'FontSize',18) 
clear PlotNames;
fflm  = 0.1:0.1:10;
PlotNames = {};
yl = [];
xl = [];
for j = 1: length(TestName)
    plot(ffFilt{j}, Alpha{j} , '.', 'Color',Cols{j}{1} ,'MarkerSize', 16);
    pname = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
    PlotNames{end +1} = pname;
    
    y = log(Alpha{j}(Alpha{j} > 0));
    x = log(ffFilt{j}(Alpha{j} > 0)); 
    
    yl = [yl;y];
    xl = [xl;x];
    
    mdl = fitlm(x,y);  
    b = mdl.Coefficients.Estimate(1);
    m = mdl.Coefficients.Estimate(2); 
    plot(fflm, exp(m*log(fflm) + b),':', 'Color',Cols{j}{1} , 'LineWidth', 2);
    pname = ['Model Slope : ', num2str(m)];
    PlotNames{end +1} = pname;
end
plot(ffmod, AB0, '-.r' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Conservative';

plot(ffmod, AB1, '-.b' , 'LineWidth', 3);
PlotNames{end +1} = 'Boltzman Losing Energy';

plot(ffmod, AEM, '--k' , 'LineWidth', 3);
PlotNames{end +1} = '2dEMM';

mdlL = fitlm(xl,yl);  
bL = mdlL.Coefficients.Estimate(1);
mL = mdlL.Coefficients.Estimate(2); 
plot(fflm, exp(mL*log(fflm) + bL),':k', 'LineWidth', 2);
PlotNames{end +1} = ['All Data Model Slope : ', num2str(mL)];

xlim([0.5 , 2.55 ]);

xlabel('Frequency');
ylabel('Attenuation Coefficient');
title(['Attenuation Coefficient', Description]);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend(PlotNames);

%Save Matrices
TpTarg = TpA;
HsTarg = HsA;
F = ffFilt;
Trans = Tran;

MatFile_NM = strcat('Data/Gen/NewTrans',int2str(conc{2}),WaveType(1:3));
save(MatFile_NM,'TpTarg','HsTarg','F','Trans');

return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return
