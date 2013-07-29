% fn_CompareAttnModels
%
% INPUTS:
%
% pers = vector of periods
% ATTN = plot (dimensional) attenuation coefficient?
% EN   = plot transmitted energy?
% DIRS = plot transmitted energy as function of angle?
% EVALS = plot eigenvalues of periodic problem?
% PRBS = [3d real geom, 3d random phases, 2d, Boltzmann steady]
% COMM = comments on (1) or off (0)

function Main_CompareAttnModels

if ~exist('COMM','var'); COMM=0; end

if ~exist('PRBS','var'); PRBS=[0,0,1,0]; end

if ~exist('col','var'); mkr='o'; end
if ~exist('col_2d','var'); col_2d='b'; end
if ~exist('col_3d','var'); col_3d='r'; end
if ~exist('col_3db','var'); col_3db='g'; end
if ~exist('col_Bol','var'); col_Bol='m'; end

if ~exist('ATTN','var'); ATTN=1*101; end
if ~exist('EN','var'); EN=0*102; end
if ~exist('DIRS','var'); DIRS=0*300; end
if ~exist('EVALS','var'); EVALS=0*401; end

if ~exist('pers','var'); pers = 0.6:0.2:2; end

if ATTN; figure(ATTN); hold on; end
if EN; figure(EN); hold on; end
if EVALS; figure(EVALS); hold on; end
if DIRS; 
 for loop_p=1:length(pers) 
  figure(DIRS+loop_p); hold on; 
  title(['period = ' num2str(pers(loop_p))])
  set(gca,'yscale','log'); 
  set(gca,'ylim',[0,1])
 end 
end

if ~exist('conc','var'); conc=.79; end

if conc==.79
 file_mkr0 = 'ela_plt_Rigid16epspt005_T';
else
 file_mkr0 = 'ela_plt_Rigid8epspt1_T';
end

%% 3d Wavetank (real geometry)
if PRBS(1)

if ~exist('Rows','var'); Rows = 5; end
if ~exist('ens0','var'); ens0 = 100; end

EnTStat = zeros(6,length(pers)); attn_3d = zeros(6,length(pers));

for loop_p=1:length(pers)
    
 file_marker = [file_mkr0 int2str(floor(pers(loop_p))) 'pt' ...
     int2str(10*(pers(loop_p)-floor(pers(loop_p))))];   

 [TraStat{loop_p},EnTStat(:,loop_p),attn_3d(:,loop_p)]=...
     fn_CombineRows(file_marker,Rows,ens0,0,COMM);
 
 M(loop_p) = length(TraStat{loop_p}(1,:));
 
end

%%% PLOTS 

%%%% Attenuation coeff 

if ATTN
 figure(ATTN)   
 plot(pers,attn_3d(1,:),[col_3d mkr],'markersize',12)
 plot(pers,attn_3d(2,:),[col_3d '*'],'markersize',10)
 plot(pers,attn_3d(3,:),[col_3d '.'],'markersize',10)
 plot(pers,attn_3d(4,:),[col_3d 'x'],'markersize',10)
 title(['full perturbed (' mkr '); periodic (*); periodic eval (.); zeroth order (x)'])
end

%%%% Transmitted energy 

if EN
 figure(EN)   
 plot(pers,EnTStat(1,:),[col_3d mkr],'markersize',12)
end

%%%% Directional spectrum 

if DIRS
 for loop_p=1:length(pers)
  figure(DIRS+loop_p) 
  plot(TraStat{loop_p}(7,:),TraStat{loop_p}(1,:),['bo'],'markersize',8)
  plot(TraStat{loop_p}(7,:),TraStat{loop_p}(5,:),['bx'],'markersize',8)
  plot(TraStat{loop_p}(7,:),TraStat{loop_p}(3,:),['rs'],'markersize',8)
  plot(TraStat{loop_p}(7,:),TraStat{loop_p}(4,:),['r+'],'markersize',8)
 end 
end

%%%% eigenvalues of periodic problem 

if EVALS
 figure(EVALS)
 for loop_p=1:length(pers)
  dum(loop_p) = max(abs(imag(TraStat{loop_p}(9,:))));
 end
 plot(pers,dum,[col_3d mkr],'markersize',12); clear dum
end

end % end PRB(1)

%% 3d Wavetank (random phases between rows)

if PRBS(2)

EnTStat = zeros(3,length(pers)); attn_3d = zeros(4,length(pers));

for loop_p=1:length(pers)
    
 file_marker = [file_mkr0 int2str(floor(pers(loop_p))) 'pt' ...
     int2str(10*(pers(loop_p)-floor(pers(loop_p))))];   

 [TraStat{loop_p},EnTStat(:,loop_p),attn_3d(:,loop_p)]=...
     fn_CombineRows(file_marker,Rows,ens0,10*pi/pers(end),COMM);
 
end

%%% PLOTS 

%%%% Attenuation coeff 

if ATTN
 figure(ATTN)   
 plot(pers,attn_3d(1,:),[col_3db mkr],'markersize',12)
 plot(pers,attn_3d(2,:),[col_3db '*'],'markersize',10)
 plot(pers,attn_3d(3,:),[col_3db '.'],'markersize',10)
 plot(pers,attn_3d(4,:),[col_3db 'x'],'markersize',10)
end

%%%% Transmitted energy 

if EN
 figure(EN)   
 plot(pers,EnTStat(1,:),[col_3db mkr],'markersize',12)
end

%%%% Directional spectrum 

if DIRS
 for loop_p=1:length(pers)
  figure(DIRS+loop_p) 
  plot(TraStat{loop_p}(3,:),TraStat{loop_p}(1,:),[col_3db mkr],'markersize',8)
  plot(TraStat{loop_p}(3,:),TraStat{loop_p}(3,:),[col_3db '*'],'markersize',8)
 end 
end

end

%% 2d attenuation (long floe limit & no long floe limit)

if PRBS(3)

if ~exist('ens1','var'); ens1 = 100; end

attn_2d = zeros(1,length(pers));

for loop_p=1:length(pers)
 attn_2d(loop_p) = Main_2dAttn('Oceanide','freq',1/pers(loop_p),ens1,conc,1,COMM);
end

%%% PLOTS 

%%%% Attenuation coeff 

if ATTN
 figure(ATTN)
 plot(pers,attn_2d,'bo','markersize',12)
end

attn_2d = zeros(1,length(pers));

for loop_p=1:length(pers)
 attn_2d(loop_p) = Main_2dAttn('Oceanide','freq',1/pers(loop_p),ens1,conc,0,COMM);
end

%%% PLOTS 

%%%% Attenuation coeff 

if ATTN
 figure(ATTN)
 plot(pers,attn_2d,'rx','markersize',12)
end

end

%% Boltzmann steady

if PRBS(4)

if ~exist('th_res','var'); th_res=50; end

I_Boltz = zeros(length(pers),th_res);

for loop_p=1:length(pers)

 [I_Boltz(loop_p,:),th_vec] = ...
     Main_Boltzmann('Oceanide','freq',1/pers(loop_p),conc,th_res,COMM,0);
 
end

if DIRS 
 for loop_p=1:length(pers)
  figure(DIRS+loop_p)
  plot(cos(th_vec),I_Boltz(loop_p,:),[col_Bol '-.'])
 end 
%  for loop_p=2:length(pers)
%   plot(cos(th_vec),(loop_p-1)+0*I_Boltz(loop_p,:),'k:')
%  end
%  set(gca,'box','on')
end

end

return
