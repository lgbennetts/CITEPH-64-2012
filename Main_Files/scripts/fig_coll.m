%% fig_coll.m
%% Author: Timothy Williams
%% Date: 20141024, 16:21:01 CEST
clear;

%% ===================================================
%% inputs to Main_Trans

input_struct.conc  = 79;
input_struct.DO_PLOT  = 1;
input_struct.DO_MODEL = 1;
input_struct.what_mod = '2d EMM';
input_struct.DO_DATA  = 1;

if 0
   %% do collisions, no drag
   input_struct.collision_inputs.use_drag = 0;

   %% wave amps
   amps = 1e-2*[1;1.5;2;4;5];
   input_struct.collision_inputs.incident_amplitudes = amps;

   %% restitution coefficients
   input_struct.collision_inputs.restitution_coefficients = .9+0*amps;
   %input_struct.collision_inputs.restitution_coefficients = 0*amps;

   %% show inputs
   disp(input_struct.collision_inputs.restitution_coefficients);
else
   %% do collisions, with drag
   input_struct.collision_inputs.use_drag = 1;

   %% wave amps
   amps = 1e-2*[1;1.5;2;4;5];
   input_struct.collision_inputs.incident_amplitudes = amps;

   %% restitution coefficients
   input_struct.collision_inputs.restitution_coefficients = .9+0*amps;
   %input_struct.collision_inputs.restitution_coefficients = 0*amps;

   %% drag coefficients
   %input_struct.collision_inputs.drag_coefficients = 1.5e-2+0*amps;
   input_struct.collision_inputs.drag_coefficients = 8e-3+0*amps;
   %input_struct.collision_inputs.drag_coefficients = 1.5e-3+0*amps;

   %% type of drag law
   %input_struct.collision_inputs.drag_law = 'linear';
   input_struct.collision_inputs.drag_law = 'quadratic';

   %% show inputs
   disp(input_struct.collision_inputs.drag_coefficients);
   disp(input_struct.collision_inputs.restitution_coefficients);
end
%% ===================================================

save('out/fig_coll.mat','input_struct');

%%main call;
Main_Trans(input_struct);

%% tidy figure;
box on;
xlim([0,2.2]);
GEN_proc_fig('Wave Period, s','Transmitted Energy');
fig_name = ['out/data_comp.eps'];
saveas(gcf,fig_name,'epsc');
