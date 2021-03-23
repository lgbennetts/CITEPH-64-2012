function [user_number,user_name]  = citeph_get_user;

user_name   = getenv('LOGNAME');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tim:
if strcmp(user_name,'timill')
   %%Tim on desktop, hexagon or nansen;
   user_number = 1;
end
if strcmp(user_name,'twilliams')
   %%Tim on laptop;
   user_number = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Luke:
if strcmp(user_name,'a1612881')
   user_number = 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jordan
if strcmp(user_name,'a1229158')
   %%Tim on desktop, hexagon or nansen;
   user_number = 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
