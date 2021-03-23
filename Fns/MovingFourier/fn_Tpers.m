function T = fn_Tpers(Tp)

if Tp==.65
 T = 10;
elseif Tp==.8
 T = 10;
elseif Tp==.95
 T = 8;
elseif Tp==1.1
 T = 8;
elseif Tp==1.25
 T = 7;
elseif Tp==1.4
 T = 5;
elseif Tp==1.55
 T = 4;
elseif Tp==1.7
 T = 4;
elseif Tp==1.85
 T = 3;
elseif Tp==2
 T = 3;
else
 T = 1;
end

return