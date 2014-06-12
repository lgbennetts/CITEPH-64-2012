% function fn_sketchbasin
%
% written for paper
%
% L Bennetts 02.01.2014 / Adelaide

function fn_sketchbasin

if ~exist('fig','var'); fig=fn_getfig; end

if ~exist('ms','var'); ms=6; end

figure(fig); hold on
set(gca,'fontsize',16)

area([-2.5,2.5],[16,16],'FaceColor',0.7*[1,1,1],'LineStyle','none')

plot(-16+[0,0],[0,16],'k','linewidth',1)
plot( 16+[0,0],[0,16],'k','linewidth',1)
plot( 16.2+[0,0],[0,16],'k','linewidth',1)
plot( 16.4+[0,0],[0,16],'k','linewidth',1)

plot([-16,16.4], 16+[0,0],'k','linewidth',1)
plot([-16,16.4], 0+[0,0],'k','linewidth',1)

[xy_lhs,xy_rhs] = citeph_sensor_spots();

plot(xy_lhs(:,1),8+xy_lhs(:,2),'k.','markersize',ms)
plot(xy_rhs(:,1),8+xy_rhs(:,2),'k.','markersize',ms)

% plot(-2.5+[0,0],8*[-1,1],'k--','linewidth',1)
% plot( 2.5+[0,0],8*[-1,1],'k--','linewidth',1)

daspect([1 1 1])

set(gca,'box','off')
set(gca,'visible','off')

return