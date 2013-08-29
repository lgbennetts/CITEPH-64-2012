function out=fn_WhatTestData(conc,TYP)

if ~exist('conc','var'); conc=79; end
if ~exist('DISP','var'); DISP=0; end
if ~exist('TYP','var'); TYP='Regular'; end

%% 79% concentration

if conc==79

if DISP; cprintf(0.6*[1,1,1],'<.79 conc>:\n'); end

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

if strcmp(TYP,'Regular')
 dum=dum_reg;
elseif strcmp(TYP,'Irregular')
 dum=dum_irr;
end

clear dum_reg dum_irr

[~,inds] = sort(dum(:,2));

dum = dum(inds,:)';

if DISP; disp([['Hs [mm] : '; 'Tp [s]  : '],num2str(dum)]); end

count=1;

while count<=length(dum(1,:))
 rep = zeros(length(dum(1,:))-1);
 for loop=count+1:length(dum(1,:))
  rep(loop)=isequal(dum(:,count),dum(:,loop));
 end
 dum(:,find(rep)) = [];
 clear rep
 count=count+1;
end

HT79 = dum;

clear dum

out=HT79;

%% 39% concentration

elseif conc==39

if DISP; cprintf(0.6*[1,1,1],'<.39 conc>:\n'); end

c39_prams = conc39_testspecs();

N = length(c39_prams);

count_r=1; count_i=1;

for loop=1:N
 if strcmp(c39_prams(loop).type,'Regular')
  dum_reg(count_r,1) = 10*c39_prams(loop).wave_height;
  dum_reg(count_r,2) = 10\c39_prams(loop).period;
  count_r=count_r+1;
 elseif strcmp(c39_prams(loop).type,'Irregular')
  dum_irr(count_r,1) = 10*c39_prams(loop).wave_height;
  dum_irr(count_r,2) = 10\c39_prams(loop).period;
  count_i=count_i+1;
 end
end

if strcmp(TYP,'Regular')
 dum=dum_reg;
elseif strcmp(TYP,'Irregular')
 dum=dum_irr;
end

clear dum_reg dum_irr

[~,inds] = sort(dum(:,2));

dum = dum(inds,:)';

if DISP; disp([['Hs [mm] : '; 'Tp [s]  : '],num2str(dum)]); end

count=1;

while count<=length(dum(1,:))
 rep = zeros(length(dum(1,:))-1);
 for loop=count+1:length(dum(1,:))
  rep(loop)=isequal(dum(:,count),dum(:,loop));
 end
 dum(:,find(rep)) = [];
 clear rep
 count=count+1;
end

HT39 = dum;

out = HT39;

end % end if conc

return