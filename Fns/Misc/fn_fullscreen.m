% function fn_fullscreen(hd)
%
% sets the figure to full screen

function fn_fullscreen(hd)

if ~exist('hd','var')
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
else
 for lp=1:length(hd)
  set(hd(lp),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
 end 
end

return