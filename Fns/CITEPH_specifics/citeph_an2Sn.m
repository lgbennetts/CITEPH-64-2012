function [Sn,ff] = citeph_an2Sn(an,df)
%% citeph_an2Sn.m
%% Author: Timothy Williams
%% Date:   20130731, 15:32:03 CEST
%%
%% converts Fourier coefficients an=[a_n]: n=0,1,...,N/2-1,-N/2,-N/2+1...,-1
%% to a spectrum Sn evaluated at ff=df/2+df*(0:N/2-1),
%% where df is the freq resolution
%%
%% an = ifft(displ), where displ is the signal;
%%  displ   = sum(an.*exp(2i*pi*t*ff0))
%%   (t \in dt*(0:N-1), dt=T/N=1/N/df);
%%  ff0     = (-N/2:N/2-1)*df;
%%
%% should have var(displ)=sum(abs(an).^2)=sum(Sn*df)~\int_0^fmax S df,
%%  where fmax=N/2*df;

N  = length(an);
ff = (0:N/2-1)'*df;

if 0%%initial attempt
   Sn(1)    = abs(an(1)^2)/df;
   jj       = 2:N/2;
   Sn(jj)   = abs( an(jj).^2 )/df + abs( an(N+1-jj).^2 )/df;

else%%way that makes the integral of S equal the signal variance exactly
   jj    = 1:N/2;
   aneg  = flipud(an(jj+N/2));%% a_n: n=-1,...,-N/2
   apos  = an(jj);            %% a_n: n=0,...N/2-1
   %%
   Sn = (abs(apos.^2)+abs(aneg.^2))/df;
   ff = ff+df/2;%%move output frequency to middle of bin instead of the left;
end
