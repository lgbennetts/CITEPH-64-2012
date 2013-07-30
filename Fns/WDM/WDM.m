% function WDM(run_num)
% WDM.m  Wavelet Directional Method 
% Driver for producing WDM directional wave spectra plots
%           from a specified array of 3 or more wave staffs.
%
% VARIABLES:
%
% np       = number of staffs/probes
% [A,R]    = angle [rad] and position [m] of the wave staffs 
%            from origin; i.e. Polar coordinates of the staffs  
% [a,r]    = angle and distance between points
% data     = array with each wave staff time series as a column.
% ns       = sampling frequency for wavelet.m
% nv       = number of voices for wavelet.m
% lf,hf    = lower/upper limits of freq window
% lp, hp   =: lp=log(lf)/log(2);lp=floor(lp);
%             hp=log(hf)/log(2);hp=ceil(hp);
% thm, ths = mean and std of angles as fn of freq and time
% kkm, kks = mean and std of wavenumber as fn of freq and time
%  
% M.Donelan / Bergen 1994.11.26; Modified by A.Babanin and M.Donelan
% 2012.05 and earlier

function WDM(run_num,TYP)

if ~exist('run_num','var'); run_num=1; end
if ~exist('DEL','var'); DEL=1; end

for run=run_num  %[062]

if run > 99
 eval(['load s13',int2str(run)])
elseif run > 9
 eval(['load s130',int2str(run)])
else
 eval(['load s1300',int2str(run)])
end

if ~exist('TYP','var'); TYP=1; end

cprintf('blue',['>>> ' description '\n'])

%if ~TYP; data=data(:,1);X=X(1);Y=Y(1); end 

np = size(data,2);
npp = (np*(np-1))/2;

% Parameters required by wavelet.m
lf=.0625;
hf=1;
nv=4;

if size(data,1) > n-1
    
%data=binavg(data,2);   % binavg;  Modify as necessary to subsample time series.

% Reduce to multiple of n pts for wavelet analysis
tlen = fix(length(data)/n); data=data(1:n*tlen,:);
% Remove mean and trend.
ws = detrend(data);

if TYP

l=0;x=[];y=[];
for j=1:np-1
 for k=(j+1):np
  l=l+1;
  x(l)=X(k)-X(j);
  y(l)=Y(k)-Y(j);
 end
end

r=abs(x+1i*y);
a=atan2(y,x);

% Determine pairs with angles ~90 or 270 [degs]
l=0;rr=[];
for j=1:npp-1
 for k=(j+1):npp
   l=l+1;
   rr(l)=a(j)-a(k);
   csj(l)=cos(a(j));
   csk(l)=cos(a(k));
   snj(l)=sin(a(j));
   snk(l)=sin(a(k));
   rk(l)=r(k);
   rj(l)=r(j);
 end
end
rr=rr*180/pi;
ii=find(rr<0);
rr(ii)=rr(ii)+360;
ii=find(rr>70 & rr<110);
jj=find(rr>250 & rr<290);
ij=sort([ii jj]);

if isempty(ij)
 cprintf('red','array does not satisfy angle criteria... exiting\n')
 return
end

clear ii jj X Y x y rr a r

end

for i1=1:n:length(ws(:,1))-n+1
 lp=log(lf)/log(2);lp=floor(lp);
 hp=log(hf)/log(2);hp=ceil(hp);

 AMP = [];
 for jj=1:np
  eval(['wx',int2str(jj),' = [];'])
  eval(['B',int2str(jj),' = [];'])
  eval(['AMP',int2str(jj),' = [];'])
  % (LB comm) pick off `blocks' size 2^12 for each probe:
  x=ws(i1:i1+n-1,jj);x=x(:);
  dum_t=tm(i1:i1+n-1);dum_t=dum_t(:);
  %
  eval(['[wx',int2str(jj),',f]=wavelet(x,lp,hp,nv,ns);'])
  eval(['B',int2str(jj),'=angle(wx',int2str(jj),');'])
  eval(['AMP',int2str(jj),'=abs(wx',int2str(jj),');'])
 end
 % (LB comm) Average the amplitude over the probes:
 AMP = AMP1;
 for jj = 2:np
  eval(['AMP = AMP + AMP',int2str(jj),';'])
  eval(['clear AMP',int2str(jj),' wx',int2str(jj)])
 end
 AMP = AMP/np;
 clear AMP1 wx1 x

 mxmf=max(find(f < ns/2));
 kkm=[];thm=[];kks=[];ths=[];
 
 if TYP
 %********************************************
 for mf=1:mxmf %%% (LB comm) loop over freqs
 %********************************************

  l=0;b=[]; %%% (LB comm) difference in phases at freq mf
  for j=1:np-1
   for k=(j+1):np
    l=l+1;
    eval(['b(:,l)=B',int2str(k),'(:,mf)-B',int2str(j),'(:,mf);'])
   end
  end
  ill=find(b(:) > pi); b(ill)=b(ill)-2*pi;
  ill=find(b(:) < -pi); b(ill)=b(ill)+2*pi;
  lb=length(b);
  l=0;AA=[];
  for j=1:npp-1
   for k=(j+1):npp
    l=l+1;
    % (LB comm) determine if staff pair satisfies angle criteria:
    if length(find(ij == l)) ==0 
     AA(:,l)=zeros(lb,1); bk(:,l)=zeros(lb,1);bj(:,l)=zeros(lb,1);
    else
     AA(:,l)=rk(l)/rj(l)*b(:,j)./b(:,k);
     bk(:,l)=b(:,k);bj(:,l)=b(:,j);
    end % end if
   end % end k
  end % end j

  th=[];ll=0;
  % (LB comm) see eqns 8-9 of Donelan et al 96
  for l=ij
   ll=ll+1;
   th(:,ll)=atan2((AA(:,l)*csk(l)-csj(l)),(snj(l)-AA(:,l)*snk(l)));
   kk(:,ll)=(snk(l)*bj(:,l)/rj(l)-snj(l)*bk(:,l)/rk(l)) ./ ...
           (csj(l)*snk(l) - csk(l)*snj(l))./cos(th(:,ll)); 
  end
  th = th+(kk<0)*pi;
  kk = abs(kk);
  thm = [thm meanang(th')'];
  ths = [ths std(th')'];
  kkm = [kkm mean(kk')'];
  kks = [kks std(kk')'];
 end %  end of 'for mf=1:mxmf'  loop (frequency bin loop)
 clear AA b bj bk ill th kk j k l
 end

 % (LB comm) save amp, waveno and angle as fns of freq/time
 if i1 == 1
  t=dum_t;
  AAp=AMP;
  ddd=round(thm*180/pi);
  kkmp=kkm;
  eval(['save ../Fns/WDM/yw',int2str(run),' AAp ddd kkmp f t'])
 else
  eval(['load ../Fns/WDM/yw',int2str(run)])
  t = [t;dum_t];
  AAp=[AAp;AMP];
  kkmp=[kkmp;kkm];
  ddd=[ddd;round(thm*180/pi)];
  eval(['save ../Fns/WDM/yw',int2str(run),' AAp ddd kkmp f t'])
 end % end if i1 ...
 clear AAp ddd kkmp f dum_t
end % end for i1, i.e. blocks of n

clear B1
for jj = 2:np
 eval(['clear B',int2str(jj)])
end

clear csj csk dr hf hp kkm kks lb ll lf mf mxmf mw A R
clear n rj rk snj snk thm ths wn t i1 ij jj lp data AMP

if ~TYP
 wavetime
else
 waveplots
 wavenums
end

else
 cprintf('m',['>>> require > ' int2str(n) 'data points\n'])
end   % end of 'if length(data) > n' loop

if DEL
 eval(['delete ../Fns/WDM/yw',int2str(run) '.mat'])
 if run > 99
  eval(['delete ../Fns/WDM/s13',int2str(run) '.mat'])
 elseif run > 9
  eval(['delete ../Fns/WDM/s130',int2str(run) '.mat'])
 else
  eval(['delete ../Fns/WDM/s1300',int2str(run) '.mat'])
 end
end % end if DEL
end % end run=run_num



return