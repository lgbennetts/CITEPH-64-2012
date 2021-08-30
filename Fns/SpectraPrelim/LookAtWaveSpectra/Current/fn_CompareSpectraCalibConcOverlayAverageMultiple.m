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
if ~exist('conc','var') conc={1,39}; end %39,79 or empty (1)
if ~exist('probes','var'); probes= {11:20,11:20} ;end;     %probes= 1:20 ;end; 
if ~exist('TestName','var'); TestName= {{'14','9'}}; end ;  
%if ~exist('TestName','var'); TestName= {{'15','10'}, {'14','9'},{'16','8'}}; end ;  

% if ~exist('conc','var') conc={39,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'10','10'}, {'9','9'},{'8','8'}}; end ;  

% if ~exist('conc','var') conc={39,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'10','10'}, {'9','9'},{'8','8'}}; end ;  


%Compare reflection to transmission
% if ~exist('conc','var') conc={39,39}; end %39,79 or empty (1)
% if ~exist('probes','var'); probes= {1:10,11:20} ;end;     %probes= 1:20 ;end; 
% if ~exist('TestName','var'); TestName= {{'10','10'}, {'9','9'},{'8','8'}}; end ;  


if ~exist('Cols','var'); Cols= {'xr', '^b','sk'}; end ;

if ~exist('binnum','var');  binnum=30; end;
if ~exist('UBfactor','var');  UBfactor=1.25; end;
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
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs = 10*c_pram.wave_height;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{1},Tp,ProbeLocXY(1),WaveType); %Slower
        TindLB = fn_Tind(conc{1},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{1},TpUB,ProbeLocXY(1),WaveType);


        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time ;

        [Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,T_window);
        period = 1./ff;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));
        
        Sn1 = mean(Sn_mat(:,jj0:jj1), 2);
        %% Fourier Transform

        Periods1{i} = period;

        SnS1{i} = Sn1;
        Times1{i} = [round(tm(SI));round(tm(EI))];



        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc{2},TestName{j}{2}, probes{2}(i));

        Tp = c_pram.period/10.0;
        Hs = 10*c_pram.wave_height;

        dt = tm(2) - tm(1);

        Tind = fn_Tind(conc{2},Tp,ProbeLocXY(1),WaveType); %Slower
        
        TpLB = Tp*LBfactor;
        TpUB = Tp*UBfactor;
        Hs = 10*c_pram.wave_height;
        
        TindLB = fn_Tind(conc{2},TpLB,ProbeLocXY(1),WaveType);
        TindUB = fn_Tind(conc{2},TpUB,ProbeLocXY(1),WaveType);
        
        t0index = find(strcmp({TindLB.description},t0description));
        t0 = TindLB(t0index).time;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = TindUB(t1index).time ;
        
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
        Times2{i} = [round(tm(SI));round(tm(EI))];

        i = i +1;
    end

%     SnS1L = [SnS1{1},SnS1{2},SnS1{3},SnS1{4},SnS1{5},SnS1{6},SnS1{7},SnS1{8},SnS1{9},SnS1{10}];
%     SnS2L = [SnS2{1},SnS2{2},SnS2{3},SnS2{4},SnS2{5},SnS2{6},SnS1{7},SnS2{8},SnS2{9},SnS2{10}];

      SnS1L = [SnS1{1}];
      SnS2L = [SnS2{1}];
% 
    InitialTimeAverage1 = mean(SnS1L,2);
    InitialTimeAverage2 = mean(SnS2L,2);
    
%     PeriodsA{j} = Periods1{1};
%     TransmissionsCoeff{j} = movmean(InitialTimeAverage2,wbinwidth) ./ movmean(InitialTimeAverage1,wbinwidth);

    StepSize = 0.02;
    BinMids = 0.6:StepSize:1.2;
    BinEdges(1) = -1;
    BinEdges(2:length(BinMids)+1) = BinMids - StepSize/2;
    BinEdges(end + 1) = BinMids(end) + StepSize/2;
    BinEdges(end+1) = 20000;

    [N,edges,bin] = histcounts(Periods1{1},BinEdges);
    
%     BinMidsTrunc = BinMids(1 + N>0);

%     [N,edges,bin] = histcounts(Periods1{1},binnum);
    ITA1 = accumarray(bin(:),InitialTimeAverage1,[],@sum);
    ITA2 = accumarray(bin(:),InitialTimeAverage2,[],@sum);

    TransmissionsCoeff{j} = ITA2(2: 1+ length(BinMids))./ ITA1(2: 1+ length(BinMids));
    
    PeriodsA{j} =  BinMids;
    
    t0str = num2str(max(cellfun(@min,Times1(1:size(SnS1,1)))));
    t1str = num2str(min(cellfun(@max,Times1(1:size(SnS1,1)))));
    Lab =['Tp = ', num2str(Tp), '(s)  Hs = ',num2str(Hs),'(mm)','  ( ',t0str,' [s])',' || ', '( ',t1str,' [s])'] ;

    plot(PeriodsA{j},TransmissionsCoeff{j},Cols{j},'DisplayName',Lab, 'MarkerSize',5);
    j = j +1;
end

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

plot(model_pers,Boltzman0,'--r','DisplayName','Boltzmann Remove Scattering','LineWidth',2);
plot(model_pers,Boltzman1,'--b','DisplayName','Boltzmann','LineWidth',2);

plot(model_pers,TwoDEMM,'-.k','DisplayName','2D EMM','LineWidth',2);

title(['Transmission Coefficients - Ratio of Binned Averages of Spectral Energy, number of bins : ' num2str(binnum)]);
xlim([0.4,2]);
ylim([0,1]);
legend();


return




function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return

