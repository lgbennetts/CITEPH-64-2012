% This program takes output of WDM and summarizes in wavenumber polar space.
% and stores the wavenumber-direction spectrum ("ksp") in the folder
% "kspect".

dr=pi/180;
% set maximum wavenumber: usually to fix(2*pi/diameter of array)
maxwn=6;
kfac=120/maxwn;
figure(10*run);clf;
%
eval(['ex2=exist([''yw'',int2str(run),''.mat'']);'])
if ex2==2
run
eval(['load yw',int2str(run)])

stdh=mean(std(detrend(ws)));
varh=stdh*stdh;
wsstats = stats(ws);

clear ws data

[m,n]=size(AAp);
runl=m;
K=zeros(360,121);
ddd=mod360(ddd);
ii=find(ddd == 0); ddd(ii)=ddd(ii)+360;

%nv = Number of voices in Wavelet transformation.
fnv=f/f(1);
fnvf=find(fnv>1.9999 & fnv < 2.0001);
nv=fnvf-1;
df=f*log(2)/nv;
waven=[1:kfac*maxwn]/kfac;dk=1/kfac;wavenn=[1:kfac*maxwn+1]/kfac;
k=round(kkmp*kfac);clear kkmp;
ik=find(k > maxwn*kfac);
k(ik)=(maxwn*kfac+1)*ones(size(ik));
for kk = 1:maxwn*kfac+1
 run
 kk
 ik=[];ik=find(k==kk);
  if ~isempty(ik);
   for di = 1:360
    id=[];id=find(ddd(ik)==di);
    K(di,kk)=K(di,kk) + sum(AAp(ik(id)).^2)/nv*1.03565/runl;
   end
  end
end
K=(K./(ones(360,1)*wavenn))*kfac*180/pi;
corr=varh/(sum(sum(K).*wavenn)/kfac*pi/180);

Hs=4*sqrt(sum(sum(K).*wavenn)/kfac*pi/180);

fnsave=['kspect\ksp_',int2str(run)];
eval(['save ',fnsave,' K waven run kfac runl Hs corr wsstats'])

subplot(221)
loglog(waven,sum(K(:,1:120)).*waven*pi/180,'.-'); grid on
xlabel('wavenumber [m^{-1}]')
ylabel('spectral density [m^3]')
subplot(222)
contour(1:360,waven(1:120),K(:,1:120)'.*(waven(1:120)'.^0*ones(1,360)));grid
title('Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
subplot(223)
contour(1:360,waven(1:120),K(:,1:120)'.*(waven(1:120)'.^3*ones(1,360)));grid
title('k^3 \times Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
subplot(224)
contour(1:360,waven(1:120),K(:,1:120)'.*(waven(1:120)'.^4*ones(1,360)));grid
title('k^4 \times Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
pause(1)

clear K waven run kfac runl Hs corr AAp ddd kk df di ik k 
else
end
