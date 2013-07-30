% This program takes output of WDM and summarizes in frequency polar space,
% and stores the frequency-direction spectrum ("fsp") in the folder
% "fspect".

jjj = run;

eval(['ex=exist([''../Fns/WDM/yw'',int2str(jjj),''.mat'']);'])

if ex==2
 eval(['load ../Fns/WDM/yw',int2str(jjj)])
 jjj;  
 stdh=mean(std(detrend(ws)));
 varh=stdh*stdh;
 wsstats = stats(ws);

 AA=AAp; clear AAp

 fnv=f/f(1);
 fnvf=find(fnv>1.9999 & fnv < 2.0001);
 nv=fnvf-1;
 ddd=mod360(ddd);
 ii=find(ddd == 0); ddd(ii)=ddd(ii)+360;
 [m,n]=size(ddd);
 runl=m;
 df=f*log(2)/nv;
 
 for j=1:n %%% (LB COMM) loop over freqs
  j;
  for k=1:360 %%% (LB COMM) loop over degrees
   kk=k;
   %if k==361;kk=400;end
   ii=find(ddd(:,j)==kk);
   E(k,j)=sum(AA(ii,j).^2)/nv*1.03565/runl; %%% (LB COMM) energy=E(freq,ang) 
  end
 end
 E=(E./(ones(360,1)*(df.*f)))*180/pi;
 corr=varh/(sum(sum(E).*df.*f)*pi/180);
 Hs=4*sqrt(sum(sum(E).*df.*f)*pi/180);

 if ~DEL
  fnsave=['../WDM/fspect/fsp_',int2str(jjj)];
  eval(['save ',fnsave,' E f Hs corr runl run wsstats']) 
 end

 figure(run);clf;
 subplot(221)
 loglog(f,sum(E(:,:)).*f*pi/180,'.-'); grid on
 xlabel('frequency [Hz]')
 ylabel('spectral density [m^2/Hz]')
 subplot(222)
 contour(1:360,f(1:n),E'.*(f(1:n)'.^0*ones(1,360)));grid
 title('Energy spectrum')
 ylabel('frequency [Hz]')
 xlabel('direction [degrees]')
 subplot(223)
 contour(1:360,f(1:n),E'.*(f(1:n)'.^4*ones(1,360)));grid
 title('f^4 \times Energy spectrum')
 ylabel('frequency [Hz]')
 xlabel('direction [degrees]')
 subplot(224)
 contour(1:360,f(1:n),E'.*(f(1:n)'.^5*ones(1,360)));grid
 title('f^5 \times Energy spectrum')
 ylabel('frequency [Hz]')
 xlabel('direction [degrees]')
 pause(1)
 %eval(['print -deps fig',int2str(jjj)]) 
 clear AA ddd df ii j k m kk E f ji

else

 cprintf('red','>>> waveplots needs a file\n')
 
end

return
