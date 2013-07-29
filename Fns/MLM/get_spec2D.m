%
clear all
close all
%
% f1='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TESTS\8630.txt';
% f2='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TESTS\8631.txt';
% f3='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TESTS\8632.txt';
% f4='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TESTS\8633.txt';
%
f1='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TEST_LowPass\lp8630.txt';
f2='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TEST_LowPass\lp8631.txt';
f3='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TEST_LowPass\lp8632.txt';
f4='D:\ALESSANDRO\EXP_TRONDHEIM\DATA\TEST_LowPass\lp8633.txt';
%
ns=512;  %length of individual timeseries !!!!!!!!!!!!!!!
fs=80;
%
%
S=0;
%
%
c=0;
%
disp(f1);
s=load(f1);
%
htest=4*std(s(8000:end,22:29));
pk=find(htest<=htest(end)/10);
ppk=isempty(pk);
%
S1=0;
if ppk==1
    c=c+1;
    [m,n]=size(s);
    W=[];
    W(:,1)=s(8001:floor(end/2)*2,22)-wave_lowpass(s(8001:floor(end/2)*2,22),80,10);
    W(:,2)=s(8001:floor(end/2)*2,23)-wave_lowpass(s(8001:floor(end/2)*2,23),80,10);
    W(:,3)=s(8001:floor(end/2)*2,24)-wave_lowpass(s(8001:floor(end/2)*2,24),80,10);
    W(:,4)=s(8001:floor(end/2)*2,25)-wave_lowpass(s(8001:floor(end/2)*2,25),80,10);
    W(:,5)=s(8001:floor(end/2)*2,26)-wave_lowpass(s(8001:floor(end/2)*2,26),80,10);
    W(:,6)=s(8001:floor(end/2)*2,27)-wave_lowpass(s(8001:floor(end/2)*2,27),80,10);
    W(:,7)=s(8001:floor(end/2)*2,28)-wave_lowpass(s(8001:floor(end/2)*2,28),80,10);
    W(:,8)=s(8001:floor(end/2)*2,29)-wave_lowpass(s(8001:floor(end/2)*2,29),80,10);
    E=spec2d_exp(W,ns,0);
    S1=S1+E.E;
end
if ppk==0;disp('    SKIP...     ');end;
%
disp(f2);
s=load(f2);
%
htest=4*std(s(8000:end,22:29));
pk=find(htest<=htest(end)/10);
ppk=isempty(pk);
%
if ppk==1
    c=c+1;
    [m,n]=size(s);
    W=[];
    W(:,1)=s(8001:floor(end/2)*2,22)-wave_lowpass(s(8001:floor(end/2)*2,22),80,10);
    W(:,2)=s(8001:floor(end/2)*2,23)-wave_lowpass(s(8001:floor(end/2)*2,23),80,10);
    W(:,3)=s(8001:floor(end/2)*2,24)-wave_lowpass(s(8001:floor(end/2)*2,24),80,10);
    W(:,4)=s(8001:floor(end/2)*2,25)-wave_lowpass(s(8001:floor(end/2)*2,25),80,10);
    W(:,5)=s(8001:floor(end/2)*2,26)-wave_lowpass(s(8001:floor(end/2)*2,26),80,10);
    W(:,6)=s(8001:floor(end/2)*2,27)-wave_lowpass(s(8001:floor(end/2)*2,27),80,10);
    W(:,7)=s(8001:floor(end/2)*2,28)-wave_lowpass(s(8001:floor(end/2)*2,28),80,10);
    W(:,8)=s(8001:floor(end/2)*2,29)-wave_lowpass(s(8001:floor(end/2)*2,29),80,10);
    E=spec2d_exp(W,ns,0);
    S1=S1+E.E;
end
if ppk==0;disp('    SKIP...     ');end;
%
disp(f3);
s=load(f3);
%
htest=4*std(s(8000:end,22:29));
pk=find(htest<=htest(end)/10);
ppk=isempty(pk);
%
if ppk==1
    c=c+1;
    [m,n]=size(s);
    W=[];
    W(:,1)=s(8001:floor(end/2)*2,22)-wave_lowpass(s(8001:floor(end/2)*2,22),80,10);
    W(:,2)=s(8001:floor(end/2)*2,23)-wave_lowpass(s(8001:floor(end/2)*2,23),80,10);
    W(:,3)=s(8001:floor(end/2)*2,24)-wave_lowpass(s(8001:floor(end/2)*2,24),80,10);
    W(:,4)=s(8001:floor(end/2)*2,25)-wave_lowpass(s(8001:floor(end/2)*2,25),80,10);
    W(:,5)=s(8001:floor(end/2)*2,26)-wave_lowpass(s(8001:floor(end/2)*2,26),80,10);
    W(:,6)=s(8001:floor(end/2)*2,27)-wave_lowpass(s(8001:floor(end/2)*2,27),80,10);
    W(:,7)=s(8001:floor(end/2)*2,28)-wave_lowpass(s(8001:floor(end/2)*2,28),80,10);
    W(:,8)=s(8001:floor(end/2)*2,29)-wave_lowpass(s(8001:floor(end/2)*2,29),80,10);
    E=spec2d_exp(W,ns,0);
    S1=S1+E.E;
end
if ppk==0;disp('    SKIP...     ');end;
%
disp(f4);
s=load(f4);
%
htest=4*std(s(8000:end,22:29));
pk=find(htest<=htest(end)/10);
ppk=isempty(pk);
%
if ppk==1
    c=c+1;
    [m,n]=size(s);
    W=[];
    W(:,1)=s(8001:floor(end/2)*2,22)-wave_lowpass(s(8001:floor(end/2)*2,22),80,10);
    W(:,2)=s(8001:floor(end/2)*2,23)-wave_lowpass(s(8001:floor(end/2)*2,23),80,10);
    W(:,3)=s(8001:floor(end/2)*2,24)-wave_lowpass(s(8001:floor(end/2)*2,24),80,10);
    W(:,4)=s(8001:floor(end/2)*2,25)-wave_lowpass(s(8001:floor(end/2)*2,25),80,10);
    W(:,5)=s(8001:floor(end/2)*2,26)-wave_lowpass(s(8001:floor(end/2)*2,26),80,10);
    W(:,6)=s(8001:floor(end/2)*2,27)-wave_lowpass(s(8001:floor(end/2)*2,27),80,10);
    W(:,7)=s(8001:floor(end/2)*2,28)-wave_lowpass(s(8001:floor(end/2)*2,28),80,10);
    W(:,8)=s(8001:floor(end/2)*2,29)-wave_lowpass(s(8001:floor(end/2)*2,29),80,10);
    E=spec2d_exp(W,ns,0);
    S1=S1+E.E;
end
if ppk==0;disp('    SKIP...     ');end;
%
E.E=S1./c;
%
filename=(['spec_2D_N24_G6_MLM_512freq_2.mat']);
save(filename,'E','-v6');
