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

function [Energy] = fn_ExtractEnergy(conc,TestName,probes,Tp,TpsRegion,TestConc,t0description,t1description, TWindow)

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
        
        [Sn_mat,tm_vec,ff] = fn_JustFourierMovingWindow(tm,disp,Tp,TWindow);
        Periods = 1./ff;

        [~,jj0]=min(abs(tm_vec-t0));
        [~,jj1]=min(abs(tm_vec-t1));

        Sn = mean(Sn_mat(:,jj0:jj1), 2);
        Time = [tm_vec(jj0), tm_vec(jj1)];
        

        
%         [~,jjLB]=min(abs(Periods-TpsRegion(1)));
%         [~,jjUB]=min(abs(Periods-TpsRegion(2)));   
        Energys(i) = sum(Sn(and(Periods >  TpsRegion(1), Periods <  TpsRegion(2) )));

% 
%         figure();
%         plot(Periods,Sn);
%         hold on;
%         plot([Periods(jjLB),Periods(jjLB)],[0,max(Sn)]);
%         plot([Periods(jjUB),Periods(jjUB)],[0,max(Sn)]);

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


