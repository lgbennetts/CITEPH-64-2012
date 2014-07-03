% function out=fn_WhatTestData(conc,TYP)
%
% DESCRIPTION: produces the periods, heights, wavelengths (non-dim wrt floe diam)
%              and steepnesses for the expts
%
% INPUTS:
%
% conc = for the attenuation tests either 39 or 79
%        for the single floe tests use conc=1
%
% OUTPUTS:
%
% out  = array of [Heights; periods]
% 
% FLAGS:
%
% DISP = diplay on (1) or 0 (off, default)
% TYP  = 'Regular' (default) or 'Irregular'
%
% L Bennetts July 2013 / La Seyne

function out=fn_WhatTestData(conc,TYP,DISP)

if ~exist('conc','var'); conc=79; end
if ~exist('TYP','var'); TYP='Regular'; end
if ~exist('DISP','var'); DISP=0; end

%% single floe tests

if conc==1
 
 if DISP; cprintf(0.6*[1,1,1],'<single floe>:\n'); end
 
 prams = singlefloe_testspecs();
 
 N = length(prams);
 
 count_r=1; count_i=1;
 
 for loop=1:N
  if strcmp(prams(loop).type,'Regular')
   dum_reg(count_r,1) = 10*prams(loop).wave_height;
   dum_reg(count_r,2) = 10\prams(loop).period;
   count_r=count_r+1;
  elseif strcmp(prams(loop).type,'Irregular')
   dum_irr(count_r,1) = 10*prams(loop).wave_height;
   dum_irr(count_r,2) = 10\prams(loop).period;
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
 
 HT = dum;
 
 clear dum
 
 out=HT;


%% 79% concentration

elseif conc==79
 
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
 
 if strfind(TYP,'Regular')
  dum=dum_reg;
 elseif strfind(TYP,'Irregular')
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
  if strfind(c39_prams(loop).type,'Regular')
   dum_reg(count_r,1) = 10*c39_prams(loop).wave_height;
   dum_reg(count_r,2) = 10\c39_prams(loop).period;
   count_r=count_r+1;
  elseif strfind(c39_prams(loop).type,'Irregular')
   dum_irr(count_r,1) = 10*c39_prams(loop).wave_height;
   dum_irr(count_r,2) = 10\c39_prams(loop).period;
   count_i=count_i+1;
  end
 end
 
 if strfind(TYP,'Regular')
  dum=dum_reg;
  if strfind(TYP,'calib'); dum([3:4,8:9],:)=[]; end
 elseif strfind(TYP,'Irregular')
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

Param = ParamDef_Oceanide;

for loop=1:size(out,2)
 Forcing = Force_def(Param.g(1), Param.bed, 'freq', 1/out(2,loop));
 out(3,loop) = Forcing.lam0/.99;
 out(4,loop) = out(3,loop)\out(1,loop)/10;
end

clear Forcing

return