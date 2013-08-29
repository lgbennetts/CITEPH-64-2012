% function Main_SingleDisk
%
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam0 = value of forcing fortyp (vector)
%
% L Bennetts Aug 2013 / Adelaide

function Main_SingleDisk

if ~exist('th_vec','var'); th_vec=linspace(-pi,pi,201); end

if ~exist('TEST','var'); TEST='Oceanide'; end
if ~exist('PLT','var');  PLT=1; end
if ~exist('COMM','var');  COMM=0; end
col='c';

%% Preliminaries

if strcmp(TEST,'Oceanide')
    
 if ~exist('RIGID','var'); RIGID=6; end   

 if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.495,33e-3]; end
 if ~exist('Param','var'); Param = ParamDef3d_Oceanide(GeomDisk,RIGID); 
    Param = ModParam_def(Param,1,10,0,0); end

 if ~exist('fortyp','var'); fortyp='freq'; end
 if ~exist('lam0','var'); lam0=1./(.25:0.05:1.85); end

 if ~exist('bed','var'); bed=3.1; end
 
else
 cprintf('green', 'Code another test!!! \n')
 return
end

%% Calculate RAOs

RAO = zeros(3,length(lam0));

for loop_lam=1:length(lam0)

 out = fn_ElasticDisk(fortyp, lam0(loop_lam), Param, GeomDisk, bed, ...
  th_vec, RIGID, 1, COMM, 0);

 for loop_out=1:length(out)
  if strcmp(out(loop_out).name,'RAO-heave')
   RAO(1,loop_lam)=out(loop_out).value;
  elseif strcmp(out(loop_out).name,'RAO-pitch')
   RAO(2,loop_lam)=out(loop_out).value;
  elseif strcmp(out(loop_out).name,'RAO-surge')
   RAO(3,loop_lam)=out(loop_out).value;
  end
 end % end loop_out
 
end % end for lam0

%% Plot

if PLT
 figure(PLT);
 if isempty(get(gcf,'children'))
  str={'RAO heave', 'RAO pitch', 'RAO surge'};
  for loop=1:size(RAO,1)
   h(loop)=subplot(1,size(RAO,1),loop); hold on
   xlabel(h(loop),'period [s]','fontsize',14)
   ylabel(h(loop),str{loop},'fontsize',14)
   set(h(loop),'fontsize',14)
  end
 else
  h=flipud(get(gcf,'children'));
 end % end if isempty(get(gcf,'children'))

 for loop=1:size(RAO,1)
  plot(h(loop),1./lam0,RAO(loop,:),col)
 end
  
end % end PLT


return