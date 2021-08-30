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

function fn_TransmissionNewAvg(conc,TestName,probes)

% %Compare calibration to disk
% if ~exist('conc','var') conc={0,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  
% if ~exist('TestName','var'); TestName= {{'14','9'},{'15','10'},{'16','8'}}; end ;  

if ~exist('conc','var') conc={0,79}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','19'}}; end ; 
if ~exist('TestName','var'); TestName= {{'14','19'}, {'15','20'},{'16','21'},{'17','22'}}; end ; 

% if ~exist('conc','var') conc={0,0}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'14','14'}}; end ;  
 


%Compare reflection to transmission
% if ~exist('conc','var') conc={39,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'10','10'}, {'9','9'},{'8','8'}}; end ;  


% if ~exist('Cols','var'); Cols= {'xr', '^b','sk','dg'}; end ;
if ~exist('Cols','var'); Cols= {{'#ff0000','#680303'},{'#0bff01 ','#058000'},{'#0487f9','#00427c'},{'#9701ff','#470179'}}; end ;



if ~exist('Vert_Modes','var'); Vert_Modes=1e2; end
if ~exist('model_pers','var'); model_pers=0.3:0.02:2; end
if ~exist('DO_FDSP','var');  DO_FDSP=0; end

if ~exist('PerNum','var');  PerNum=7; end

if ~exist('TpersStr','var');  TpersStr='Tpers=fn_Tpers(Tp,WaveType);'; end
if ~exist('OutputStr','var');  OutputStr = 'TargetSpectra'; end

%{'14','9'}
%{'15','10'}
%{'16','8'}
close all;

for j = 1: length(TestName)
    
    [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{1},TestName{j}{1}, 1);
    Tp = c_pram.period/10.0;
    Hs = c_pram.wave_height/100.0;
    
    TpS = linspace(0.75*Tp,Tp,PerNum);
    TpS = [TpS(1:end-1),linspace(Tp,1.5*Tp,PerNum)];
    
    PerWholeCalib = [];
    SpectraWholeCalib = [];
    PerWholeConc = [];
    SpectraWholeConc = [];
    for ji = 1: length(TpS)
        OutputStr = 'TargetSpectra';
%         outCalib = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{1},TestName{j}{1},probes{1},OutputStr);
        outCalib = fn_SpectraAverageTgT(TpS(ji),conc{1},TestName{j}{1},WaveType,probes{1},OutputStr);
        PerFilt = outCalib(1).value;
        SpectraFilt = outCalib(2).value;
        
        PerWholeCalib = [PerFilt;PerWholeCalib];
        SpectraWholeCalib = [SpectraFilt;SpectraWholeCalib];
        
%         outConc = fn_SpectraAroundTgT(TpC,TpS(ji-1),TpS(ji),conc{2},TestName{j}{2},probes{2},OutputStr);
        outConc = fn_SpectraAverageTgT(TpS(ji),conc{2},TestName{j}{2},WaveType,probes{2},OutputStr);
        PerFilt = outConc(1).value;
        SpectraFilt = outConc(2).value;
        
        PerWholeConc = [PerFilt;PerWholeConc];
        SpectraWholeConc = [SpectraFilt;SpectraWholeConc];        
    end

    ffWholeConc = 1./PerWholeConc;
    
    wJS = 2*pi/Tp;
    ffJS = 0:0.01:10;
    JSspec{j}  = jonswap(2*pi*ffJS,'wp',wJS,'Hs',Hs);
    
    %Figure
     figure();
     plot(PerWholeCalib,SpectraWholeCalib,'-b', 'DisplayName','Calibration');
     hold on;
     plot(PerWholeConc,SpectraWholeConc,'-r', 'DisplayName','Experiment');
     plot(1./ffJS,JSspec{j},'--k' , 'DisplayName','JonSwap');
     title(['Tp = ', num2str(Tp), ' Hs = ', num2str(Hs)]);
     xlabel('T(s)')
     ylabel('Spectra')
    
   
    %experimental
    TpA{j} = Tp;
    HsA{j} = Hs;
    ffA{j} = ffWholeConc;
    ffpA{j} =ffWholeConc*Tp;
    
   
    ATWCalib{j} = SpectraWholeCalib;
    ATWConc{j} = SpectraWholeConc;
    
%     ATWCalibmax{j} = max(SnS1L,[],2);
%     ATWCalibmin{j} = min(SnS1L,[],2);
%     ATWConcmax{j} = max(SnS2L,[],2);
%     ATWConcmin{j} = min(SnS2L,[],2);
    
    
    fpL = 0.75;
    fpH = 1.5;
    AmmTol = 10^-3;
    STol = AmmTol^2/2;
%     STol = 0;
    
    %above 1mm accuracy
    Cond1 = ATWConc{j} > STol;
    Cond2 = and(ffpA{j} > fpL, ffpA{j} < fpH);
   
    CondInd = and(Cond1,Cond2);
    
    ffFilt{j} = ffA{j}(CondInd);
    ffpFilt{j} = ffpA{j}(CondInd);
    ATWCalibFilt{j} = ATWCalib{j}(CondInd);    
    ATWConcFilt{j} = ATWConc{j}(CondInd);
    
%     ATWCalibFiltmax{j} = ATWCalibmax{j}(CondInd);    
%     ATWConcFiltmax{j} = ATWConcmax{j}(CondInd);
%     
%     ATWCalibFiltmin{j} = ATWCalibmin{j}(CondInd);    
%     ATWConcFiltmin{j} = ATWConcmin{j}(CondInd);
    
    Tran{j} = (ATWConcFilt{j}./ATWCalibFilt{j});
%     TranMin{j} = (ATWConcFiltmin{j}./ATWCalibFilt{j});
%     TranMax{j} = (ATWConcFiltmax{j}./ATWCalibFilt{j});
    
    Alpha{j} = -log(Tran{j})./5;
    clear wlc;
    for i = 1: length(ffFilt{j})
        ki = dispersion_free_surface((2*pi*ffFilt{j}(i))^2/9.81,0,3.1);
        wlc(i,1) = 2*pi/imag(ki);
    end
    
    wlExp{j} = wlc;
    
    j = j +1;
end

%Save Matrices
TpTarg = TpA;
HsTarg = HsA;
F = ffFilt;
WL =  wlExp;
Trans = Tran;
Ain = ATWCalibFilt;

MatFile_NM = strcat('Data/Gen/TransEnergy',int2str(conc{2}),WaveType(1:3));
save(MatFile_NM,'TpTarg','HsTarg','F','WL','Trans','Ain');

% 
% Description = ['Concentration ' num2str(conc{2}) '% '];
% 
% %frequency plot range - smooth
% ffJSsmoothnn = 0:0.01:100;

%Plot both 1.4 tests!, then comparison of Experiment, Calibration and
%Target
%Then transmission with filter and attenuation with filter

%good evidence of linearity in regions, is both 1.4 tests
% figure();
% set(gca,'FontSize',18) 
% hold on;
% PlotNames = {};
% for j = 1:length(TpA)
%     pTp = 1.4;
%     wJS = 2*pi/pTp;
%     JSSmooth  = jonswap(2*pi*ffJSsmoothnn,'wp',wJS,'Hs',HsA{1});
%     JSSmoothNorm = JSSmooth / (HsA{1}/4)^2;
%     ffJSSmooth = ffJSsmoothnn*pTp;
%     
%     if TpA{j} == 1.4
%         m0Calib = trapz(ffA{j},ATWCalib{j});
% 
%         plot(ffpA{j}, ATWCalib{j}./m0Calib , '-', 'Color',Cols{j}{1} , 'LineWidth', 3);
%         pnameCalib = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%         PlotNames{end +1} = pnameCalib;
% 
% %         Upper = ATWCalibmax{j}(2:end)./m0Calib;
% %         Lower = ATWCalibmin{j}(2:end)./m0Calib;
% %         ffrem0 = ffpA{j}(2:end);
% %         ShadeX = [ffrem0; flipud(ffrem0)];
% %         inBetween = [Upper ; flipud(Lower)];
% %         patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
% %         alpha(0.3)
% %         PlotNames{end +1} = 'Max and Min Region (Among Probes)';
% % 
% %         plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %         PlotNames{end +1} = 'Max and Min Boundary';
% %         plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %         PlotNames{end +1} = '';
%     end
%         
%     
% end
% plot(ffJSSmooth, JSSmoothNorm, '--k' , 'LineWidth', 3);
% PlotNames{end +1} = 'JonSwap [m]';

% xlabel('Frequency / Fp');
% ylabel('Spectra / m0 ');
% title(['Calibration Spectra Comparisons', Description]);
% legend(PlotNames);
% xlim([0 , 3 ]);
% 
% %Copy Figure, with log scales
% a1 = gca;
% f2 = figure;
% a2 = copyobj(a1,f2) ;
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlim([10^-1 , 10^2 ]);
% ylim([10^-8 , 10 ]);
% legend(PlotNames);

%Plot of Calibration vs Experiment Comparisons

%Spectra Comparisons - both Calibration and Experimental

% for j = 1: length(TestName)
%     figure();
%     set(gca,'FontSize',18) 
%     hold on;
%     PlotNames = {};
%     
%     m0Calib = trapz(ffA{j},ATWCalib{j});
%     
%     plot(ffpA{j}, ATWCalib{j}./m0Calib , '-', 'Color',Cols{j}{1} , 'LineWidth', 3);
%     pnameCalib = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%     PlotNames{end +1} = pnameCalib;
%     
% %     Upper = ATWCalibmax{j}(2:end)./m0Calib;
% %     Lower = ATWCalibmin{j}(2:end)./m0Calib;
% %     ffrem0 = ffpA{j}(2:end);
% %     ShadeX = [ffrem0; flipud(ffrem0)];
% %     inBetween = [Upper ; flipud(Lower)];
% %     patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
% %     alpha(0.5)
% %     PlotNames{end +1} = 'Max and Min Region (Among Probes)';
% %     
% %     plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %     PlotNames{end +1} = 'Max and Min Boundary';
% %     plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %     PlotNames{end +1} = '';
% 
%     plot(ffpA{j}, ATWConc{j}./m0Calib , '-', 'Color',Cols{j}{2} , 'LineWidth', 3);
%     pnameConc = ['Conc ', num2str(conc{2}),' | Probes : ( ', num2str(probes{2}(1)),' , ', num2str(probes{2}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%     PlotNames{end +1} = pnameConc;
%     
% %     Upper = ATWConcmax{j}(2:end)./m0Calib;
% %     Lower = ATWConcmin{j}(2:end)./m0Calib;
% %     ffrem0 = ffpA{j}(2:end);
% %     ShadeX = [ffrem0; flipud(ffrem0)];
% %     inBetween = [Upper ; flipud(Lower)];
% %     patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
% %     alpha(0.5)
% %     PlotNames{end +1} = 'Max and Min Region (Among Probes)';
% %     
% %     plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %     PlotNames{end +1} = 'Max and Min Boundary';
% %     plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %     PlotNames{end +1} = '';
%     
%     wJS = 2*pi/TpA{j};
%     JSSmooth  = jonswap(2*pi*ffJSsmoothnn,'wp',wJS,'Hs',HsA{j});
%     JSSmoothNorm = JSSmooth / (HsA{j}/4)^2;
%     ffJSSmooth = ffJSsmoothnn*TpA{j};
%     
%     plot(ffJSSmooth, JSSmoothNorm, '--k' , 'LineWidth', 3);
%     PlotNames{end +1} = 'JonSwap [m]';
% 
%     xlabel('Frequency / Fp');
%     ylabel('Spectra / m0 ');
%     title(['Calibration Spectra Comparisons', Description]);
%     legend(PlotNames);
%     xlim([0 , 3 ]);
% 
%     a1 = gca;
%     f2 = figure;
%     a2 = copyobj(a1,f2) ;
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
%     xlim([10^-1 , 10^2 ]);
%     ylim([10^-8 , 10 ]);
%     legend(PlotNames);
% end



% 
% % Models
% if conc{2}==39
% dum_c = 100*pi*(0.495^2)/2;
% elseif conc{2}==79
% dum_c = 100*pi*(0.495^2);
% end
% 
% %Predictions
% Boltzman0 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',0,Vert_Modes,DO_FDSP,0);
% Boltzman0 = (Boltzman0.value).^2;
% Boltzman0 = Boltzman0.';
% Boltzman0 = [0;Boltzman0;1];
% 
% Boltzman1 = Main_AttnModels(model_pers(2:end-1),dum_c,'Boltzmann steady',1,Vert_Modes,DO_FDSP,0);
% Boltzman1 = (Boltzman1.value).^2;
% Boltzman1 = Boltzman1.';
% Boltzman1 = [0;Boltzman1;1];
% 
% %2Demm prediction
% TwoDEMM = Main_AttnModels(model_pers(2:end-1),dum_c,'2d EMM',0,Vert_Modes,DO_FDSP,0);
% TwoDEMM = (TwoDEMM.value).^2;
% TwoDEMM = TwoDEMM.';
% TwoDEMM = [0;TwoDEMM;1];
% 
% ffmod = 1.0./ model_pers.';
% ffmod = ffmod(2:end-1);
% TEM = TwoDEMM(2:end-1);
% TB0 = Boltzman0(2:end-1);
% TB1 = Boltzman1(2:end-1);
% AEM = -(log(TEM)./5);
% AB1 = -(log(TB1)./5);
% AB0 = -(log(TB0)./5);
% 
% for i = 1: length(ffmod)
%     ki = dispersion_free_surface((2*pi*ffmod(i))^2/9.81,0,3.1);
%     wlMod(i) = 2*pi/imag(ki);
% end
% 
% %Transmission PLots
% %Spectra Comparisons - both Calibration and Experimental
% figure();
% hold on;
%  set(gca,'FontSize',18) 
% clear PlotNames;
% PlotNames = {};
% for j = 1: length(TestName)
%     plot(ffFilt{j}, Tran{j} , '.', 'Color',Cols{j}{1} ,'MarkerSize', 16);
%     pname = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%     PlotNames{end +1} = pname;
%     
% %     Upper = TranMax{j};
% %     Lower = TranMin{j};
% %     ffrem0 = ffFilt{j};
% %     ShadeX = [ffrem0; flipud(ffrem0)];
% %     inBetween = [Upper ; flipud(Lower)];
% %     patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
% %     alpha(0.3)
% %     PlotNames{end +1} = 'Max and Min Region (Among Probes)';
% %     
% %     plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %     PlotNames{end +1} = 'Max and Min Boundary';
% %     plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
% %     PlotNames{end +1} = '';
%     
% end
% plot(ffmod, TB0, '-.r' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Conservative';
% 
% plot(ffmod, TB1, '-.b' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Losing Energy';
% 
% plot(ffmod, TEM, '--k' , 'LineWidth', 3);
% PlotNames{end +1} = '2dEMM';
% 
% xlabel('Frequency');
% ylabel('Transmission Coefficient ');
% title(['Transmission Coefficient ', Description]);
% legend(PlotNames);
% xlim([0 , 5 ]);
% ylim([0 , 1 ]);


%Transmission PLots - WaveLength
%Spectra Comparisons - both Calibration and Experimental
% figure();
% hold on;
%  set(gca,'FontSize',18) 
% clear PlotNames;
% PlotNames = {};
% for j = 1: length(TestName)
%     plot(wlExp{j}, Tran{j} , '.', 'Color',Cols{j}{1} ,'MarkerSize', 16);
%     pname = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%     PlotNames{end +1} = pname;
    
%     Upper = TranMax{j};
%     Lower = TranMin{j};
%     ffrem0 = wlExp{j};
%     ShadeX = [ffrem0; flipud(ffrem0)];
%     inBetween = [Upper ; flipud(Lower)];
%     patch(ShadeX, inBetween, 'r', 'LineStyle','none', 'FaceColor',Cols{j}{1})
%     alpha(0.3)
%     PlotNames{end +1} = 'Max and Min Region (Among Probes)';
%     
%     plot(ffrem0, Upper, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
%     PlotNames{end +1} = 'Max and Min Boundary';
%     plot(ffrem0, Lower, ':', 'Color',Cols{j}{1}  , 'LineWidth', 1);
%     PlotNames{end +1} = '';
    
% end
% plot(wlMod, TB0, '-.r' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Conservative';
% 
% plot(wlMod, TB1, '-.b' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Losing Energy';
% 
% plot(wlMod, TEM, '--k' , 'LineWidth', 3);
% PlotNames{end +1} = '2dEMM';
% 
% xlabel('Wavelength');
% ylabel('Transmission Coefficient ');
% title(['Transmission Coefficient ', Description]);
% legend(PlotNames);
% xlim([0 , 5 ]);
% ylim([0 , 1 ]);


% %Attenuation Coefficient
% %What if we just fit one?
% figure();
% hold on;
% set(gca,'FontSize',18) 
% clear PlotNames;
% fflm  = 0.1:0.1:10;
% PlotNames = {};
% yl = [];
% xl = [];
% for j = 1: length(TestName)
%     plot(ffFilt{j}, Alpha{j} , '.', 'Color',Cols{j}{1} ,'MarkerSize', 16);
%     pname = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%     PlotNames{end +1} = pname;
%     
%     y = log(Alpha{j}(Alpha{j} > 0));
%     x = log(ffFilt{j}(Alpha{j} > 0)); 
%     
%     yl = [yl;y];
%     xl = [xl;x];
%     
%     mdl = fitlm(x,y);  
%     b = mdl.Coefficients.Estimate(1);
%     m = mdl.Coefficients.Estimate(2); 
%     plot(fflm, exp(m*log(fflm) + b),':', 'Color',Cols{j}{1} , 'LineWidth', 2);
%     pname = ['Model Slope : ', num2str(m)];
%     PlotNames{end +1} = pname;
% end
% plot(ffmod, AB0, '-.r' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Conservative';
% 
% plot(ffmod, AB1, '-.b' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Losing Energy';
% 
% plot(ffmod, AEM, '--k' , 'LineWidth', 3);
% PlotNames{end +1} = '2dEMM';
% 
% mdlL = fitlm(xl,yl);  
% bL = mdlL.Coefficients.Estimate(1);
% mL = mdlL.Coefficients.Estimate(2); 
% plot(fflm, exp(mL*log(fflm) + bL),':k', 'LineWidth', 2);
% PlotNames{end +1} = ['All Data Model Slope : ', num2str(mL)];
% 
% xlim([0.5 , 2.55 ]);
% 
% xlabel('Frequency');
% ylabel('Attenuation Coefficient');
% title(['Attenuation Coefficient', Description]);
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend(PlotNames);


%Attenuation Coefficient
%What if we just fit one?
% figure();
% hold on;
% set(gca,'FontSize',18) 
% clear PlotNames;
% fflm  = 0.1:0.1:10;
% for i = 1: length(fflm)
%     ki = dispersion_free_surface((2*pi*fflm(i))^2/9.81,0,3.1);
%     wllm(i) = 2*pi/imag(ki);
% end
% 
% 
% PlotNames = {};
% yl = [];
% xl = [];
% for j = 1: length(TestName)
%     plot(wlExp{j}, Alpha{j} , '.', 'Color',Cols{j}{1} ,'MarkerSize', 16);
%     pname = ['Conc ', num2str(conc{1}),' | Probes : ( ', num2str(probes{1}(1)),' , ', num2str(probes{1}(end)) ' ) | Hs ', num2str(HsA{j}),' [m] | Tp ', num2str(TpA{j}),' [s]'];
%     PlotNames{end +1} = pname;
%     
%     y = log(Alpha{j}(Alpha{j} > 0));
%     x = log(wlExp{j}(Alpha{j} > 0)); 
%     
%     yl = [yl;y];
%     xl = [xl;x];
%     
%     mdl = fitlm(x,y);  
%     b = mdl.Coefficients.Estimate(1);
%     m = mdl.Coefficients.Estimate(2); 
%     plot(wllm, exp(m*log(wllm) + b),':', 'Color',Cols{j}{1} , 'LineWidth', 2);
%     pname = ['Model Slope : ', num2str(m)];
%     PlotNames{end +1} = pname;
% end
% plot(wlMod, AB0, '-.r' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Conservative';
% 
% plot(wlMod, AB1, '-.b' , 'LineWidth', 3);
% PlotNames{end +1} = 'Boltzman Losing Energy';
% 
% plot(wlMod, AEM, '--k' , 'LineWidth', 3);
% PlotNames{end +1} = '2dEMM';
% 
% mdlL = fitlm(xl,yl);  
% bL = mdlL.Coefficients.Estimate(1);
% mL = mdlL.Coefficients.Estimate(2); 
% plot(wllm, exp(mL*log(wllm) + bL),':k', 'LineWidth', 2);
% PlotNames{end +1} = ['All Data Model Slope : ', num2str(mL)];
% 
% xlim([0.5 , 2.55 ]);
% 
% xlabel('Wavenumber');
% ylabel('Attenuation Coefficient');
% title(['Attenuation Coefficient', Description]);
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend(PlotNames);


return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

