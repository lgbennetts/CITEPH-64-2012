% function fn_whichprobe(side,probe,fig)
%
% function to identify probe in array
%
% L Bennetts 31.07.2013 / St Cyr, France
%
% INPUTS:
%
% side  = 'LHS' or 'RHS'
% probe = 1->10 (can be multiple)
% fig   = figure number (999 default) 

function fn_whichprobe(side,probe,fig)

if ~exist('side','var'); side='RHS'; end
if ~exist('probe','var'); probe=1:10; end
if ~exist('fig','var'); fig=999; end

if ~exist('ms','var'); ms=12; end

figure(fig)
set(gca,'fontsize',16)

[xy_lhs,xy_rhs,lab_lhs,lab_rhs] = citeph_sensor_spots();
if strcmp(side,'LHS')
 plot(xy_lhs(:,1),xy_lhs(:,2),'.','markersize',ms)
elseif strcmp(side,'RHS')
 plot(xy_rhs(:,1),xy_rhs(:,2),'.','markersize',ms)
end

daspect([1 1 1])
hold on

if strcmp(side,'LHS')
 for loop_p=1:length(probe)
  plot(xy_lhs(probe(loop_p),1),xy_lhs(probe(loop_p),2),'m.','markersize',ms)
  title([side ' probe ' int2str(probe(loop_p)) ': ' lab_lhs(loop_p)])
  pause
  dum=get(gca,'children');
  delete(dum(1))
 end
elseif strcmp(side,'RHS')
 for loop_p=1:length(probe)
  plot(xy_rhs(probe(loop_p),1),xy_rhs(probe(loop_p),2),'m.','markersize',ms)
  title([side ' probe ' int2str(probe(loop_p)) ': ' lab_rhs(loop_p)])
  pause
  dum=get(gca,'children');
  delete(dum(1))
 end
end

return

