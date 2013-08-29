function Main_RAO

if ~exist('file_pre','var'); file_pre = 'Temp_data/s00'; end

if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_FPLT','var');  DO_FPLT=0; end
if ~exist('DO_PLOT','var');  DO_PLOT=1; end

HT=fn_WhatTestData(1,'Regular',0);

ht_inds=1:length(HT);

probes=1;
rbms={'surge'};

%% INC WAVE

for loop_ht=ht_inds
 for probe=probes
  Main_SingleFloeData(HT(2,loop_ht),HT(1,loop_ht),loop_ht,probe,...
   'WP',DO_FPLT,'none',file_pre)
 end
end

%% RIGID BODY MODES

for loop_ht=ht_inds
 for loop_rbm=1:length(rbms)
  probe=rbms{loop_rbm};
  Main_SingleFloeData(HT(2,loop_ht),HT(1,loop_ht),loop_ht,probe,...
   'RBM',DO_FPLT,'none',file_pre)
 end
end

%% RESPONSE AMPLITUDE OPERATORS

%%% Calculate incident amplitude as average over individual probes

inc_amp=zeros(length(ht_inds),1);
count_ht=1;
for loop_ht=ht_inds
 file_nm=fn_get_filenm(loop_ht,file_pre);
 eval(['load ' file_nm ' outputs;'])
 for probe=probes
  out=outputs{probe};
  inc_amp(count_ht)=inc_amp(count_ht)+out.value;
  clear out
 end
 clear outputs
 count_ht=count_ht+1;
end
inc_amp = inc_amp/length(probes);

%%% RAOs

count_ht=1;
for loop_ht=ht_inds
 file_nm=fn_get_filenm(loop_ht,file_pre);
 eval(['load ' file_nm ' outputs;']) 
 for loop_rbm=1:length(rbms)
  probe=rbms{loop_rbm};
  out=outputs{length(probes)+loop_rbm};
  RAO(count_ht,loop_rbm)=out.value;
  clear out
  if or(or(strcmp('surge',probe),strcmp('sway',probe)),strcmp('heave',probe))
   RAO(count_ht,loop_rbm)=RAO(count_ht,loop_rbm)/inc_amp(count_ht);
  end
 end
 clear outputs
 count_ht=count_ht+1;
end

%% DELETE

if DEL
 for loop_ht=ht_inds
  file_nm=fn_get_filenm(loop_ht,file_pre);
  eval(['delete ' file_nm '.mat'])
 end
end

%% PLOT

fig=200;

for loop_rbm=1:length(rbms)
 figure(fig+loop_rbm)
 plot(HT(2,ht_inds),RAO(:,loop_rbm))
end
    
return

%%% SUBFUNCTIONS %%%

function file_nm=fn_get_filenm(run,file_pre)

if ~exist('file_pre','var'); file_pre = 'Temp_data/s13'; end

if run > 99
  file_nm = [file_pre int2str(run)];
 elseif run > 9
  file_nm = [file_pre '0' int2str(run)];
 else
  file_nm = [file_pre '00' int2str(run)];
end

return 