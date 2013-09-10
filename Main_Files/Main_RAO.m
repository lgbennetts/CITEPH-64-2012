% function Main_RAO
%
% DESCRIPTION: Compare model & data RAOs for a single disk
%
% L Bennetts Aug 2013 / Adelaide

function Main_RAO

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL

if ~exist('DO_FDSP','var');  DO_FDSP=0; end
if ~exist('DO_PLOT','var');  DO_PLOT=1; end
if ~exist('DO_DISP','var');  DO_DISP=0; end

if DO_PLOT
 if ~exist('col','var'); 
  col=' ''r.'' , ''markersize'' , 12'; 
  col_model=' ''c'' ';
 end
end

%% DATA

if ~exist('DO_DATA','var'); DO_DATA=0; end

if ~exist('file_pre','var'); file_pre = 'Temp_data/s00'; end

if ~exist('T_pers','var'); T_pers=15; end

if ~exist('data_out','var') 
 data_out.name='amp-harmo-steady-1st';
 data_out.tint='t0=t0+4*Tp; t1=t0+10*Tp;'; 
end 

if ~exist('DEL','var');      DEL=1; end
if ~exist('DO_FPLT','var');  DO_FPLT=0; end 
%if ~exist('DO_FPLT','var');  DO_FPLT='Aspec-signal'; end

HT=fn_WhatTestData(1,'Regular',0); ht_inds=1:length(HT); %

probes=[1,2,3,8,9]; %1:10;
rbms={'heave','surge+sway','roll+pitch'};

%% MODEL 

if ~exist('DO_MODEL','var'); DO_MODEL=1; end

if ~exist('Vert_Modes','var'); Vert_Modes=10; end
if DO_DISP; model_pers=HT(2,ht_inds); end
if ~exist('model_pers','var'); model_pers=0.3:0.05:1.85; end

count_rbm=1;
for loop_rbm=1:length(rbms)
 if strfind(rbms{loop_rbm},'heave')
  model_outs{count_rbm} = 'RAO-heave';
 elseif strfind(rbms{loop_rbm},'surge')
  model_outs{count_rbm} = 'RAO-surge';
 elseif strfind(rbms{loop_rbm},'pitch')
  model_outs{count_rbm} = 'ANG-pitch';
 end
 count_rbm=count_rbm+1;
end
clear count_rbm
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% NUMERICAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_MODEL
 RAO_model = Main_SingleDiskModel('freq',1./model_pers,...
  Vert_Modes,model_outs,0,DO_FDSP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_DATA

%% INC WAVE

for loop_ht=ht_inds
 for probe=probes
  Main_SingleDiskData(HT(2,loop_ht),HT(1,loop_ht),T_pers,data_out,...
   loop_ht,probe,'WP',DO_FPLT,'none',DO_FDSP,file_pre)
  if DO_FPLT
   pause
   close all
  end
 end
end

%% RIGID BODY MODES

for loop_ht=ht_inds
 for loop_rbm=1:length(rbms)
  if strfind(rbms{loop_rbm},'+')
   ind = strfind(rbms{loop_rbm},'+');
   probe=rbms{loop_rbm}(1:ind-1);
   Main_SingleFloeData(HT(2,loop_ht),HT(1,loop_ht),T_pers,data_out,...
    loop_ht,probe,'RBM',DO_FPLT,'none',DO_FDSP,file_pre)
   probe=rbms{loop_rbm}(ind+1:end);
   Main_SingleFloeData(HT(2,loop_ht),HT(1,loop_ht),T_pers,data_out,...
    loop_ht,probe,'RBM',DO_FPLT,'none',DO_FDSP,file_pre)
  else
   probe=rbms{loop_rbm};
   Main_SingleFloeData(HT(2,loop_ht),HT(1,loop_ht),T_pers,data_out,...
    loop_ht,probe,'RBM',DO_FPLT,'none',DO_FDSP,file_pre)
  end % END IF '+'
  if DO_FPLT
   pause
   close all
  end
 end
end

%% RESPONSE AMPLITUDE OPERATORS (RAOs)

%%% Calculate incident amplitude as average over individual probes

if isempty(probes)
 inc_amp=ones(length(ht_inds),1);
else
 inc_amp=zeros(length(ht_inds),1);
 count_ht=1;
 for loop_ht=ht_inds
  file_nm=fn_get_filenm(loop_ht,file_pre);
  eval(['load ' file_nm ' outputs;'])
  count_pb=1;
  for probe=probes
   out=outputs{count_pb};
   inc_amp(count_ht)=inc_amp(count_ht)+out.value;
   clear out
   count_pb=count_pb+1;
  end
  clear outputs count_pb
  count_ht=count_ht+1;
 end
 inc_amp = inc_amp/length(probes);
end % END IF ISEMPTY(PROBES)

%%% RAOs

count_ht=1;
for loop_ht=ht_inds
 file_nm=fn_get_filenm(loop_ht,file_pre);
 eval(['load ' file_nm ' outputs;'])
 count_rbm=1;
 for loop_rbm=1:length(rbms)
  probe=rbms{loop_rbm};
  if strfind(probe,'+')
   out=outputs{length(probes)+count_rbm};
   dum_RAO=out.value;
   count_rbm=count_rbm+1;
   out=outputs{length(probes)+count_rbm};
   RAO(count_ht,loop_rbm)=abs(out.value + 1i*dum_RAO);
   count_rbm=count_rbm+1; clear dum_RAO
  else
   probe=rbms{loop_rbm};
   out=outputs{length(probes)+count_rbm};
   RAO(count_ht,loop_rbm)=out.value;
   count_rbm=count_rbm+1;
  end % END IF '+'
  if 0
   if or(or(~isempty(strfind(probe,'surge')),...
    ~isempty(strfind(probe,'sway'))),...
    ~isempty(strfind(probe,'heave')))
    RAO(count_ht,loop_rbm)=RAO(count_ht,loop_rbm)/inc_amp(count_ht);
   end
  else
   RAO(count_ht,loop_rbm)=RAO(count_ht,loop_rbm)/inc_amp(count_ht);
  end
  clear out
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

end % IF DO_DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%  PLOT & DISPLAY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT

if DO_PLOT
 fig=200;
 for loop_rbm=1:length(rbms)
  figure(fig+loop_rbm)
  if DO_DATA
   eval(['plot(HT(2,ht_inds),RAO(:,loop_rbm),' col ')'])
  end 
  xlabel('period [s]','fontsize',14)
  ylabel('RAO','fontsize',14)
  title(rbms{loop_rbm},'fontsize',14)
  if DO_MODEL
   hold on
   eval(['plot(model_pers,RAO_model(loop_rbm,:),' col_model ')'])
   hold off
  end % END IF DO_MODEL
 end % END LOOP_RBM
end % END DO_PLOT

%% DISPLAY

if DO_DISP
 cprintf(0.4*[1,1,1],['Periods:\n'])
 disp(HT(2,ht_inds))
 for loop_rbm=1:length(rbms)
  if and(~exist('k0','var'),...
    or(or(or(~isempty(strfind(rbms{loop_rbm},'surge')),...
     ~isempty(strfind(rbms{loop_rbm},'sway'))),...
     ~isempty(strfind(rbms{loop_rbm},'pitch'))),...
     ~isempty(strfind(rbms{loop_rbm},'roll'))))
   bed=3.1; k0=zeros(1,length(ht_inds));
   count_ht=1;
   for loop_ht=ht_inds
    sig = ((2*pi/HT(2,loop_ht))^2)/9.81;
    k0(count_ht) = CalcRealRoot_PWC([bed, sig], ...
     'FS_DispRel_PWC', 'UppLimReal_FS_PWC', 1e-16);
    count_ht=count_ht+1;
   end
   clear count_ht sig
  end % end calc k0
 end
 
 for loop_rbm=1:length(rbms)
  count_ht=1;
  for loop_ht=ht_inds
   if or(~isempty(strfind(rbms{loop_rbm},'surge')),...
     ~isempty(strfind(rbms{loop_rbm},'sway')))
   disp_vec(:,count_ht)=[RAO(count_ht,loop_rbm);...
    coth(k0(count_ht)*bed)];
  elseif or(~isempty(strfind(rbms{loop_rbm},'pitch')),...
    ~isempty(strfind(rbms{loop_rbm},'roll')))
   disp_vec(:,count_ht)=[RAO(count_ht,loop_rbm)*inc_amp(count_ht);...
    tan(pi*RAO(count_ht,loop_rbm)*inc_amp(count_ht)/180)/k0(count_ht);...
    RAO(count_ht,loop_rbm)];
  else
   disp_vec(:,count_ht)=[RAO(count_ht,loop_rbm)];
  end
  count_ht=count_ht+1;
  end % end loop_ht
  if or(~isempty(strfind(rbms{loop_rbm},'surge')),...
     ~isempty(strfind(rbms{loop_rbm},'sway')))
   cprintf(0.4*[1,1,1],[rbms{loop_rbm} ' RAO, eccentricity:\n'])
  elseif or(~isempty(strfind(rbms{loop_rbm},'pitch')),...
    ~isempty(strfind(rbms{loop_rbm},'roll')))
   cprintf(0.4*[1,1,1],[rbms{loop_rbm} ' angle (degs), tan(ang)/k, angle/amp:\n'])
  else
   cprintf(0.4*[1,1,1],[rbms{loop_rbm} ' RAO:\n'])
  end
  disp(disp_vec)
  if DO_MODEL
   cprintf(0.4*[1,1,1],['Model ' model_outs{loop_rbm} ':\n'])
   disp(RAO_model(loop_rbm,:))
  end % END IF DO_MODEL
  clear disp_vec
 end % for loop_rbm
end % end DO_DISP

return

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

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