function basedir = citeph_user_specifics()
%% citeph_user_specifics.m
%% Author: Timothy Williams
%% Date:   20130731, 09:46:06 CEST
%% basedir contains results_preliminary & calibration_waves directories with data in them

user_name   = getenv('LOGNAME');
if strcmp(user_name,'timill')
   basedir  = '/work/timill/CITEPH-data/results_preliminary/';
   cprintf('blue','Tim: your basedir needs modification\n'); 
elseif strcmp(user_name,'a1612881')
 dum = {'/Volumes/Scratch/Data/','/Volumes/My_Passport/'};
 eval(['ex=exist(dum{1});'])
 count=1;
 while ex~=7
  count=count+1;
  if count<length(dum)
   cprintf('green',['>>> ' dum ' not found\n'])
   s = input('Enter new base dir or <return> to exit', 's');
   if isempty(s)
    return
   else
    eval(['ex=exist(s);'])
   end
  else
   eval(['ex=exist(dum{count});'])
  end
 end
 basedir  = [dum{count} 'CITEPH-64-2012/'];
end

return