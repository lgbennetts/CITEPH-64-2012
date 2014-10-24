%% fig_coll.m
%% Author: Timothy Williams
%% Date: 20141024, 16:21:01 CEST
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conc  = 79;
name_str = ' ''conc'' ';
val_str  = '   conc   ';

DO_PLOT  = 1;
name_str = [name_str,';', ' ''DO_PLOT'' '];
val_str  = [val_str, ';', '   DO_PLOT   '];

DO_MODEL = 1;
name_str = [name_str,';', ' ''DO_MODEL'' '];
val_str  = [val_str, ';', '   DO_MODEL   '];

what_mod = '2d EMM';
name_str = [name_str,';', ' ''what_mod'' '];
val_str  = [val_str, ';', '   what_mod   '];

DO_DATA  = 1;
name_str = [name_str,';', ' ''DO_DATA'' '];
val_str  = [val_str, ';', '   DO_DATA   '];

if 0
   wave_amp = 5e-2;
else
   wave_amp = 1e-2*[1;1.5;2;4;5];
end
rest_coeff        = .05;
collision_inputs  = [wave_amp,rest_coeff+0*wave_amp];
name_str          = [name_str,';', ' ''collision_inputs'' '];
val_str           = [val_str, ';', '   collision_inputs   '];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['input_struct = struct(''Name'',{',...
       name_str,'}, ''Value'', {',...
       val_str,'});'] );

if 1
   for j=1:length(input_struct)
      disp(input_struct(j));
   end
   %return;
end

%%main call;
Main_Trans(input_struct);

%% tidy figure;
box on;
xlim([0,2.2]);
GEN_proc_fig('Wave Period, s','Transmitted Energy');
fig_name = ['out/data_comp_rc',num2str(rest_coeff),'.eps'];
saveas(gcf,fig_name,'epsc');
