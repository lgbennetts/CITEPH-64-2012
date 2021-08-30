function T = fn_Tpers(Tp,WaveType)

if strcmp(WaveType,'Regular')
    if Tp==.65
     T = 10;
    elseif Tp==.8
     T = 10;
    elseif Tp==.95
     T = 16;
    elseif Tp==1.1
     T = 8;
    elseif Tp==1.25
     T = 6;
    elseif Tp==1.4
     T = 6;
    elseif Tp==1.55
     T = 4;
    elseif Tp==1.7
     T = 4;
    elseif Tp==1.85
     T = 4;
    elseif Tp==2
     T = 4;
    else
     T = 1;
    end
    T = 10;
% %     T = 20;
else
%     if Tp==.8
%      T = 15;
%     elseif Tp==1.4
%      T = 10;
%     elseif Tp==2
%      T = 5;
%     else
%      T = 1;
%     end 
    T = 40;
end
return