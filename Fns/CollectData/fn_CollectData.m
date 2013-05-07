% function fn_CollectData(file_marker,fig,col)
%
% displays transmitted energy and 1st mode reflected energy

function fn_CollectData(file_marker,fig,col)

if ~exist('file_marker','var'); file_marker='0'; end
if ~exist('col','var'); col={'b','k','r'}; end

path_root = '../../../../../Documents/MatLab/Data/3d_Wavetank/Tests/';

% pwd
% ls(path_root)

file_name = ['Main_Channel_freq_',file_marker,'.mat'];

load([path_root,file_name],'lam_vec','R_vec','T_vec','v_vecs','w','reson_mkr')

warning off
load([path_root,file_name],'clockout','duration','v_info')
warning on

if exist('clockout','var'); 
 display(['Created: ' int2str(clockout(3)),'/', int2str(clockout(2)),'/',int2str(clockout(1)),' ' ,...
  int2str(clockout(4)),':', int2str(clockout(5)),':',int2str(round(clockout(6)))] ); 
end
if exist('duration','var'); 
 display(['Lasted: ' int2str(round(duration(1))),' hours ', ...
     int2str(round(duration(2))),' mins ',int2str(round(duration(3))) ' secs'] ); 
end
if exist('v_info','var'); disp(v_info); end

%% - ENERGY CHECK & Tranmitted energy- %%
for loop_lam=1:length(lam_vec)
 disp_error = []; mkr = 1;
 real_inds = find(real(v_vecs{loop_lam})~=0); Y0_Dim = length(real_inds);
 Vert_Dim = size(R_vec{loop_lam},1)/length(v_vecs{loop_lam});
 epsm = ones(1,length(real_inds)); epsm(1) = 2; epsm = epsm/2;
 for loop_Y=1:Vert_Dim:Y0_Dim*Vert_Dim
  EngErr = epsm(loop_Y)*v_vecs{loop_lam}(real_inds(loop_Y)) ...
      - sum(epsm.*v_vecs{loop_lam}(real_inds).*...
       ( abs(R_vec{loop_lam}(real_inds,real_inds(loop_Y)).').^2 ...
       + abs(T_vec{loop_lam}(real_inds,real_inds(loop_Y)).').^2 ) );

  if abs(EngErr)>1e-5
   if mkr 
    disp_error{1} = ['Energy error ' 'wk/pi = ' num2str(w*lam_vec(loop_lam)/pi)]; 
    disp_error{2} = '; Modes: '; disp_error{3} = ' -> '; mkr=0;
   end
   disp_error{2} = [disp_error{2} int2str(loop_Y) ', '];
   disp_error{3} = [disp_error{3} num2str(abs(EngErr))  ', '];
  end   
 end
 
 if ~isempty(disp_error); disp([disp_error{1},disp_error{2},disp_error{3}]); end
 
 %%% Transmitted energy
 En(1,loop_lam) = sum( transpose(epsm.*v_vecs{loop_lam}(real_inds)) .* ...
     ( squeeze(abs(T_vec{loop_lam}(real_inds,real_inds(1))).^2) ) );
end

%% Reflection & transmission coeffs

for loop_lam=1:length(lam_vec)
 R(loop_lam) = R_vec{loop_lam}(1,1);
 %T(:,loop_lam) = T_vec{loop_lam}(:,1);
 M(loop_lam) = length(R_vec{loop_lam}(:,1));
end

%% Plot

if ~exist('fig','var')
 figure
else
 if isempty(fig)
  figure
 else
  figure(fig);
 end
end

for loop=1:2
 h(loop) = subplot(2,1,loop); hold on
end

if ~exist('w','var'); w=16; end

plot(h(1),w*lam_vec/pi,abs(R).^2,col{1})
plot(h(1),w*lam_vec/pi,reson_mkr.*abs(R).^2,col{1},'linestyle','none','marker','o')
ylabel(h(1),'|R|^2')
% plot(h(2),w*lam_vec/pi,M,col{2}); set(h(2),'ylim',[0,max(M)])
% ylabel(h(2),'No modes')
plot(h(2),w*lam_vec/pi,En,col{3})
plot(h(2),w*lam_vec/pi,reson_mkr.*En,col{3},'linestyle','none','marker','o')
ylabel(h(2),'Trans energy')
xlabel(h(2),'w*k/\pi')

return