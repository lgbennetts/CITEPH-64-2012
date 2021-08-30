% function fn_JustFourierMovingWindow
%
% DESCRIPTION: Generate the transmissions matrices in Main/Data
%
% INPUTS:
% tm - time series times in seconds
% data - data at times (such as displacement in meters)
% Tpers - string containing function call, to be evaluated to give Tpers.
%
% General:
%
%
%
% Jordan Pitt - Adelaide - 2021 - based on old version of Main_Trans by
% Luke Bennets - 2013. 

function [Energy] = fn_ExtractEnergy(conc,TestName,probes,Tp,TpsRegion,TestConc,t0description,t1description)

    i = 1;
    
    for probei= probes
        [ProbeLocXY,tm,disp,c_pram,WaveType,Success] = fn_FindAndReadProbe(conc,TestName, probei);
        Tind = fn_Tind(TestConc,Tp,ProbeLocXY(1),WaveType);
        
        TindLB  = fn_Tind(TestConc,TpsRegion(1),ProbeLocXY(1),WaveType);
        TindUB  = fn_Tind(TestConc,TpsRegion(2),ProbeLocXY(1),WaveType);
        
        t0index = find(strcmp({TindLB.description},t0description));
        t0 = Tind(t0index).time;
        t1index = find(strcmp({TindUB.description},t1description));
        t1 = Tind(t1index).time;  
        
        t0 = max(1,t0);
        t1 = min(tm(end),t1);
        th = (t1 + t0) /2;
%         TindLB = fn_Tind(conc,TpsRegion(1),ProbeLocXY(1));
%         TindUB = fn_Tind(conc,TpsRegion(2),ProbeLocXY(1));
        dt = tm(2) - tm(1);

%         t0 = TindLB(7).time;
%         t1 = TindUB(8).time;

        dtlength = (t1 - t0) / dt;

        timewindowlength = 2^floor(log2(dtlength));
        
%         SI = floor((t0 -tm(1)) / dt);
%         EI = SI + timewindowlength-1;
        
        MI = floor((th -tm(1)) / dt);
        SI = MI - timewindowlength/2;
        EI = MI + timewindowlength/2-1;
        
        
        [Sn,mt,ff] = fn_JustFourierTransform(tm(SI:EI),disp(SI:EI));
        periods = 1./ff;
        
%         [~,jj] = min( abs(periods - Tp));
%         
%         a=Sn(1:jj-2);
%         b=Sn(2:jj-1);
%         c=Sn(3:jj);
%         TpCILB  = find((b<a & b<c)| (a ==0 | b ==0 | c == 0),1,'last')+1;
% 
%         if isempty(TpCILB)
%             TpCILB = 1;
%         end
% 
% 
%         %To Right of TPCI in period space
%         a=Sn(jj:end-2);
%         b=Sn(jj+1:end-1);
%         c=Sn(jj+2:end);
%         TpCIRB  = find((b<a & b<c)| (a ==0 | b ==0 | c == 0),1,'first')+jj;
% 
%         if isempty(TpCIRB)
%             TpCIRB = length(Sn) ;
%         end
        
        [~,jjLB]=min(abs(periods-TpsRegion(1)));
        [~,jjUB]=min(abs(periods-TpsRegion(2)));   
        Energys(i) = sum(Sn(jjUB:jjLB));
%         Energys(i) = sum(Sn(and(periods >= TpsRegion(1), periods <= TpsRegion(2))));
%         Energys(i) = sum(Sn(TpCILB : TpCIRB));
%         if i == 1
%         figure();
%         plot(periods,Sn);
%         hold on;
%         plot([periods(TpCILB),periods(TpCILB)],[0,max(Sn)]);
%         plot([periods(TpCIRB),periods(TpCIRB)],[0,max(Sn)]);
%         end
        i = i + 1;
    end
    
    Energy = mean(Energys);
    
 
    


return







function Tind = fn_Tind(conc,Tp,X,WaveType)

if conc == 39 || conc == 79
  Tind = fn_TestTimes(1/Tp,X,'attn',WaveType);
else 
  Tind = fn_TestTimes(1/Tp,X,'calibration',WaveType);
end

 return


