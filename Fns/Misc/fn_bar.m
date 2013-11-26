function fn_bar(X,T,w,hd,col,lw)

if exist('hd','var')
 if and(~isempty(hd),hd)
  axes(hd); 
 else
  axes(gca);
 end
else
 axes(gca); 
end

if ~exist('col','var'); col = 'r'; end
if ~exist('lw','var'); lw = 1.0; end
if ~exist('tol','var'); tol = 0; end

for lp=1:length(X)
 plot(X(lp)+[0,0],T(:,lp),'color',col,'linewidth',lw)
 plot(X(lp)+[-w,w],T(1,lp)+[0,0],'color',col,'linewidth',lw)
 plot(X(lp)+[-w,w],T(2,lp)+[0,0],'color',col,'linewidth',lw)
end

return