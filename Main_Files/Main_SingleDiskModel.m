% function Main_SingleDiskModel
%
% INPUTS:
%
% fortyp     = 'freq' or 'wlength' or 'waveno'
% lam0       = value of forcing fortyp (vector)
% Vert_Modes = number of vertical modes used to solve scattering problem 
% outputs    = list of the outputs required
%              'E' 
%              'E0'
%              'E_Fourier'
%              'RAO-surge'
%              'RAO-heave'
%              'RAO-pitch'
%              'RAO-pitch/k'
%              'ANG-pitch'
% FLAGS:
%
% TEST = 'Oceanide' for expts at BGO 1st July 2013
%         ... or code something else!
% PLT  = plot the output (default=on)
% COMM = turn comments on/off (off by default)
%
% OUTPUTS:
%
% RAO = array of RAOs with rows corresponding to variable outputs
%
% L Bennetts Aug 2013 / Adelaide

function RAO = Main_SingleDiskModel(fortyp,lam0,Vert_Modes,outputs,PLT,COMM)

if ~exist('outputs','var');outputs={'RAO-surge','RAO-heave','ANG-pitch'};end

if ~exist('th_vec','var'); th_vec=linspace(-pi,pi,201); end

if ~exist('Vert_Modes','var'); Vert_Modes=10; end

if ~exist('TEST','var'); TEST='Oceanide'; end
if ~exist('PLT','var');  PLT=1; col='b'; end
if ~exist('COMM','var');  COMM=0; end

%% Preliminaries

if strcmp(TEST,'Oceanide')
    
 if ~exist('RIGID','var'); RIGID=5; end   

 if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.495,33e-3]; end
 if ~exist('Param','var'); Param = ParamDef3d_Oceanide(GeomDisk,RIGID); 
    Param = ModParam_def(Param,1,Vert_Modes,0,0); end

 if ~exist('fortyp','var'); fortyp='freq'; end
 if ~exist('lam0','var'); lam0=1./(.3:0.1:1.85); end

 if ~exist('bed','var'); bed=3.1; end
 
else
 cprintf('green', 'Code another test!!! \n')
 return
end

%% Calculate RAOs

RAO = zeros(length(outputs),length(lam0));

for loop_lam=1:length(lam0)

 out = fn_ElasticDisk(fortyp, lam0(loop_lam), Param, GeomDisk, bed, ...
  th_vec, RIGID, 1, COMM, 0);

 for loop_out=1:length(out)
  for loop_outs=1:length(outputs)
   if strcmp(out(loop_out).name,outputs{loop_outs})
    str{loop_outs}=out(loop_out).name;
    RAO(loop_outs,loop_lam)=out(loop_out).value;
   end
  end % end loop_outs
 end % end loop_out
 
end % end for lam0

%% Plot

if PLT
 figure(PLT);
 if isempty(get(gcf,'children'))  
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
 
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
  
end % end PLT


return