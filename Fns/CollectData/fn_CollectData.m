function fn_CollectData(file_marker,fig,col)

if ~exist('file_marker','var'); file_marker='0'; end
if ~exist('col','var'); col={'b','k','r'}; end

path_root = '../../../../../Documents/MatLab/Data/3d_Wavetank/Tests/';

% pwd
% ls(path_root)

file_name = ['Main_Channel_freq_',file_marker,'.mat'];

load([path_root,file_name],'lam_vec','R_vec','T_vec','v_vecs','w','reson_mkr')

%% Reflection & transmission coeffs

for loop_lam=1:length(lam_vec)
 R(loop_lam) = R_vec{loop_lam}(1,1);
 %T(:,loop_lam) = T_vec{loop_lam}(:,1);
 M(loop_lam) = length(R_vec{loop_lam}(:,1));
end

%% Transmitted energy
 
for loop_lam=1:length(lam_vec)
 real_inds = find(real(v_vecs{loop_lam})~=0); 
 epsm = ones(1,length(real_inds)); epsm(1) = 2; epsm = epsm/2;

 En(1,loop_lam) = sum( transpose(epsm.*v_vecs{loop_lam}(real_inds)) .* ...
     ( squeeze(abs(T_vec{loop_lam}(real_inds,real_inds(1))).^2) ) );
end

if ~exist('fig','var')
 figure
else
 if isempty(fig)
  figure
 else
  figure(fig);
 end
end

for loop=1:3
 h(loop) = subplot(3,1,loop); hold on
end

if ~exist('w','var'); w=16; end

plot(h(1),w*lam_vec/pi,abs(R).^2,col{1})
plot(h(1),w*lam_vec/pi,reson_mkr.*abs(R).^2,col{1},'linestyle','none','marker','o')
ylabel(h(1),'|R|^2')
plot(h(2),w*lam_vec/pi,M,col{2}); set(h(2),'ylim',[0,max(M)])
ylabel(h(2),'No modes')
plot(h(3),w*lam_vec/pi,En,col{3})
plot(h(3),w*lam_vec/pi,reson_mkr.*En,col{3},'linestyle','none','marker','o')
ylabel(h(3),'Trans energy')
xlabel(h(3),'w*k/\pi')

return