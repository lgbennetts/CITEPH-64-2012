% function fn_WhatData(flpre)
%
% display the help data or selected data (with prefix flpre) 

function fn_WhatData(flpre)

if ~exist('flpre','var'); flpre=[]; end

path_root = '../../../../../Documents/MatLab/Data/3d_Wavetank/Tests/';

cat = dir([path_root ...
    flpre '*']);

for loop=1:length(cat)
 if cat(loop).name(1)~='.';
  load([path_root,cat(loop).name],'description')
  disp(cat(loop).name)
  for loop2=1:length(description)
   disp(description(loop2))
  end
 end
end

return