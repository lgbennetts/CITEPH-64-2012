% wavetime

jjj=run;

eval('ex=exist([''../Fns/WDM/yw'',int2str(jjj),''.mat'']);')

if ex==2

 eval(['load ../Fns/WDM/yw',int2str(jjj)],'AAp f')
 
 contour(t,f,AAp.')
 xlabel('time [s]')
 ylabel('freq [Hz]')
 colorbar
 
end