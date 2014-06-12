function fn_bar(X,T,w,FT,hd,col)

if ~exist('X','var'); X=1; end
if ~exist('T','var'); T=[-1;1]; end
if ~exist('w','var'); w=0.1; end

if exist('hd','var')
 if and(~isempty(hd),hd)
  axes(hd); 
 else
  axes(gca);
 end
else
 axes(gca); 
end

if ~exist('FT','var'); FT=0; end

if ~exist('col','var'); col = ' ''color'' , 0.7*[1,1,1] '; end

for lp=1:length(X)
 eval(['plot(X(lp)+[0,0],T(:,lp),' col ')'])
 if FT
  eval(['plot(X(lp)+[-w,w],T(1,lp)+[0,0],' col ')'])
  eval(['plot(X(lp)+[-w,w],T(2,lp)+[0,0],' col ')'])
 else
  if T(1,1)<T(2,1)
   T=flipud(T);
  end
  eval(['plot(X(lp)+[-w,0],T(1,lp)+[-w,0],' col ')'])
  eval(['plot(X(lp)+[ 0,w],T(1,lp)+[0,-w],' col ')'])
  eval(['plot(X(lp)+[-w,0],T(2,lp)+[w,0],' col ')'])
  eval(['plot(X(lp)+[ 0,w],T(2,lp)+[0,w],' col ')'])
 end
end

return