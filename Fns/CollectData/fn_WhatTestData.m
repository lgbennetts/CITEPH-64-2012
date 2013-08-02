function [HT79,HT39]=fn_WhatTestData

cprintf(0.6*[1,1,1],'<.79 conc>:\n')

c79_prams = conc79_testspecs();

N = length(c79_prams);

count_r=1; count_i=1;

for loop=1:N
 if strcmp(c79_prams(loop).type,'Regular')
  dum_reg(count_r,1) = 10*c79_prams(loop).wave_height;
  dum_reg(count_r,2) = 10\c79_prams(loop).period;
  count_r=count_r+1;
 elseif strcmp(c79_prams(loop).type,'Irregular')
  dum_irr(count_r,1) = 10*c79_prams(loop).wave_height;
  dum_irr(count_r,2) = 10\c79_prams(loop).period;
  count_i=count_i+1;
 end
end

[~,inds] = sort(dum_reg(:,2));

dum_reg = dum_reg(inds,:)';

disp([['Hs [mm] : '; 'Tp [s]  : '],num2str(dum_reg)])

count=1;

while count<=length(dum_reg(1,:))
 rep = zeros(length(dum_reg(1,:))-1);
 for loop=count+1:length(dum_reg(1,:))
  rep(loop)=isequal(dum_reg(:,count),dum_reg(:,loop));
 end
 dum_reg(:,find(rep)) = [];
 clear rep
 count=count+1;
end

HT79 = dum_reg;

clear dum_reg

cprintf(0.6*[1,1,1],'<.39 conc>:\n')

c39_prams = conc39_testspecs();

N = length(c39_prams);

for loop=1:N
 if strcmp(c39_prams(loop).type,'Regular')
  dum_reg(count_r,1) = 10*c79_prams(loop).wave_height;
  dum_reg(count_r,2) = 10\c79_prams(loop).period;
  count_r=count_r+1;
 elseif strcmp(c39_prams(loop).type,'Irregular')
  dum_irr(count_r,1) = 10*c79_prams(loop).wave_height;
  dum_irr(count_r,2) = 10\c79_prams(loop).period;
  count_i=count_i+1;
 end
end

[~,inds] = sort(dum_reg(:,2));

dum_reg = dum_reg(inds,:)';

disp([['Hs [mm] : '; 'Tp [s]  : '],num2str(dum_reg)])

count=1;

while count<=length(dum_reg(1,:))
 rep = zeros(length(dum_reg(1,:))-1);
 for loop=count+1:length(dum_reg(1,:))
  rep(loop)=isequal(dum_reg(:,count),dum_reg(:,loop));
 end
 dum_reg(:,find(rep)) = [];
 clear rep
 count=count+1;
end

HT39 = dum_reg;

return