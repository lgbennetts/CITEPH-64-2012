% function fig=fn_getfig
%
% DESCRIPTION: find the next figure handle available
%
% L Bennetts Sept 2013 / Adelaide
% 
% updated Nov 2014 for MatLab 2014b

function fig=fn_getfig(num)

if ~exist('num','var'); num=1; end

figHandles = findobj('Type','figure');

if isempty(figHandles)
 fig=1:num;
 if num==1; fig=fig(1); end
else
 for lp=length(figHandles):-1:1
  dum_fh(length(figHandles)+1-lp) = get(figHandles(lp),'Number');
 end
 dum = 1:max(dum_fh)+num;
 dum(dum_fh)=[];
 fig = dum(1:num);
 if num==1; fig=fig(1); end
end

return