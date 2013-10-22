function [Sn,ff,an_cx,t_mean,xtra] = citeph_get_ft(time,displ,data_type)
%% CALL: [Sn,ff,an,t_mean,xtra] = citeph_get_ft(time,displ,data_type)
%% time,displ are input vectors
%% data_type is 1 (elevation) or 2 (acceleration - need to integrate F. series wrt time)
%% Sn is the spectrum corresponding to the positive frequencies ff
%%  (ff(1)=0);
%% an are the Fourier coefficients corresponding to [ff(1:end);-fmax;-flipud(ff(end:-1:2))]
%%  (ie no fftshift applied)
%% t_mean   = mean(time)
%% xtra  = {Hs,Tp,T_m02}, T_m02=sqrt(m0/m2);

%% high pass filter by subtracting the mean;
displ = displ-mean(displ);

N        = length(time);
dt       = time(2)-time(1);%%time resolution
T        = time(N)+dt-time(1);%%record length
t_mean   = mean(time);
%%
df    = 1/T;%%freq resolution
fs    = 1/dt;%%sample rate
fmax  = fs/2;%%Nyqvist freq (max possible freq)
%%
ff0   = (-N/2:N/2-1)'*df;
%%
an       = ifft(displ,N);

if 0%data_type==2
   %%was getting strange results from this
   om    = 2*pi*fftshift(ff0);
   an    = an./(-om.^2);
   an(1) = 0;
end

jj       = 1:N/2;
aneg     = flipud(an(jj+N/2));%%-1,...,-N/2
apos     = [an(jj);aneg(N/2)];%%0,...N/2 (a_{N/2}=a_{-N/2})
an_cx    = 2*apos;%% complex amplitudes
                  %% - it seems for a real transform aneg contains no extra information
                  %% ( we just have a_{-n}=conj(a_n) );
a0       = apos(1);

if 0
   jt = 1:10;
   %[apos(jt),aneg(jt)]
   an_cx(jt)
   pause
end

if 0%%stay on LH endpoint of interval; don't use n=\pm N/2 modes
   ff    = (0:N/2-1)'*df;
   Sn    = 0*ff;
   %%
   Sn(1)    = abs(an(1)^2)/df;
   jj       = 2:N/2-1;
   Sn(jj)   = .5*abs( an_cx(jj).^2 )/df;
elseif 1%%stay on LH endpoint of interval,but do use n=\pm N/2 modes
   ff    = (0:N/2)'*df;
   Sn    = 0*ff;
   %%
   Sn(1)    = abs(a0^2)/df;
   jj       = 2:N/2;
   Sn(jj)   = .5*abs( an_cx(jj).^2 )/df;
else%%use mid-points of intervals, and use the n=-N/2 mode
   ff    = (0:N/2-1)'*df+df/2;
   Sn    = 0*ff;
   %%
   jj = 1:N/2;
   aneg  = flipud(an(jj+N/2));%%-1,...,-N/2
   apos  = an(jj);%%0,...N/2-1
   Sn    = (abs(apos.^2)+abs(aneg.^2))/df;
end
%je       = N/2+1;
%Sn(je)   = 2*abs( an(je).^2 )/df;
%Sn = [abs(an(1))^2;abs(an(2:N2/2)).^2+abs(an(N2:-1:N2/2+2)).^2]/df;

%ff       = (0:df:fmax-df)';
periods  = 1./ff;
%%
[Smax,jmax] = max(Sn);
Tp          = periods(jmax);

m0    = sum(Sn)*df;
m2    = sum(ff.^2.*Sn)*df;
T_m02 = sqrt(m0/m2);%% NB don't need 2\pi factor since moments are wrt freqency.

Hs    = 4*sqrt(m0);
xtra  = {Hs,Tp,T_m02};
%pause

if 0
   tst_var  = [var(displ,1), sum(abs(an.^2)), sum(Sn)*df]
   pause(1);
end
