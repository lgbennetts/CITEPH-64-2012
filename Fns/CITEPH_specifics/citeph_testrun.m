%% citeph_testrun.m
%% Author: Timothy Williams
%% Date:   20130716, 15:28:12 CEST

[ID_lhs,ID_rhs,freqs,nfft_lhs,nfft_rhs,wave_TH,Hs_ind]   = citeph_IDstruct();
Nprobes  = length(ID_lhs.datatypes)

method_num        = 2;
[EP_lhs,SM_lhs]   = citeph_SM_EPstruct(freqs,nfft_lhs,method_num);
[EP_rhs,SM_rhs]   = citeph_SM_EPstruct(freqs,nfft_rhs,method_num);

options  = {'MESSAGE',2,'PLOTTYPE',2,'FILEOUT',''};

addpath ..

cc = clock;
disp(['start time: ' num2str(cc(4)) ':' num2str(cc(5))]);
warning off;
%%
t0                      = rem(now,1);
%[SM_out_lhs,EP_lhs2]    = dirspec(ID_lhs,SM_lhs,EP_lhs,options);
[SM_out_rhs,EP_rhs2]    = dirspec(ID_rhs,SM_rhs,EP_rhs,options);
dt                      = round(24*60*(rem(now,1)-t0));%%mins
disp(['time taken: ' num2str(dt) ' mins']);
wave_TH
Hs_ind

if 0
   save SMs.mat ID_lhs SM_out_lhs EP_lhs2 SM_out_rhs EP_rhs2 ID_rhs
end

if 1
   SS       = SM_out_rhs.S;
   nd       = size(SS,2);
   dth      = 360/nd;
   %%
   SS0      = real(sum(SS*dth,2));%%integrate over direction
   periods  = 1./freqs;
   [Sm,jm]  = max(SS0);
   Tp       = periods(jm)
   %%
   figure(5);
   plot(periods,SS0);
   xlim([0 5]);
   GEN_proc_fig('Period, s','Frequency spectrum, m^2s');
end
