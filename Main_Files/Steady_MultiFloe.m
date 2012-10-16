function Steady_MultiFloe

disp('%---------- START: Steady_MultiFloe -----------%')

%% Parameters

Vert_Dim=2;     % - Vertical modes 
lam0 = 1.7*pi;    % - wavelength
evs=1;          % - no. evanescent waves (horiz) to include
res = 100;      % - for the integration of the Green's fns
extra_pts = 0;  % - for irreg freqs

%% Definition of the geometry of the problem solved
 GeomPlate=[2.5 2.2, 2.0, 1e-1];      %2.40482555769577
%            6.0 1.6, 0.6, 1e-1];
%            4.0 0.6, 0.4, 1e-1];
%[4.5 2.0, 1,   9e-6 ];
%           7   2.7, 0.8, 9e-1];
%           5.8 3.7, 0.6, 9e-1]; % - centre (x,y), radius, thickness

% GeomPlate=zeros(10,4);
% GeomPlate(:,2)=2.5;
% GeomPlate(:,3)=0.9;
% GeomPlate(:,4)=1e-1;
% GeomPlate(:,1)=[0:2:18].';
       
%[glob,R_max]=global_system(GeomPlate);

[Geom_Vec, thicks, rads, cx, cy] = fn_GeomDef3d(GeomPlate);

%% Check that plates don't overlap &/or touch side walls

distij=zeros(Geom_Vec(2));
angij=zeros(Geom_Vec(2));

for loop1=1:Geom_Vec(2)
 for loop2=1:Geom_Vec(2)
  [distij(loop1,loop2),angij(loop1,loop2)]=coordinate_change(...
      [cx(loop1),cy(loop1)], [cx(loop2),cy(loop2)]);        
 end
end

for loop1=1:Geom_Vec(2)-1
 for loop2=loop1+1:Geom_Vec(2)
  %distij = sqrt((cx(loop1)-cx(loop2))^2+(cy(loop1)-cy(loop2))^2);
  sumRads = rads(loop1)+rads(loop2);
   if distij(loop1,loop2)<sumRads
     disp(['plates ',num2str(loop1),' & ',num2str(loop2),' overlapping!!!!!'])
     return
   end
   clear sumRads %distij 
 end
end

if min(GeomPlate(:,2)-GeomPlate(:,3))<0
 [mn,pl]=min(GeomPlate(:,2)-GeomPlate(:,3))
 disp(['plate ',num2str(pl),' overlaps y=0'])
end

if max(GeomPlate(:,2)+GeomPlate(:,3))>Geom_Vec(5)
 [mn,pl]=max(GeomPlate(:,2)+GeomPlate(:,3))
 disp(['plate ',num2str(pl),' overlaps y=w'])
end

%% Definition of the physical parameters 
parameter_vector = ParamDef3d_v2(Geom_Vec,thicks);

freq = FindFreq_FS([inf,inf,inf,inf,parameter_vector(1)],2*pi/lam0,Geom_Vec(3));

% T0=2; % Period
% freq = 2*pi/T0;
kappa = freq^2/parameter_vector(1,1);

%kappa = 0.3385;

%% Inititializing the frequency dependent variables

x_length = Geom_Vec(4); width = Geom_Vec(end);
x_res = 101; y_res = 51;

x_vec = linspace(0,x_length,x_res);
y_vec = linspace(0,width,y_res);

%%% CALCULATE THE MESH POINTS IN THE FREE SURF REGION %%%
% - The x vals that bound the plate covered region
minx0 = min(GeomPlate(:,1)-GeomPlate(:,3));
maxx0 = max(GeomPlate(:,1)+GeomPlate(:,3));

FS_mesh = zeros(y_res,x_res);

for loop_x=1:length(x_vec)
 if and(x_vec(loop_x)>minx0,x_vec(loop_x)<maxx0)
  for loop_y=1:y_res
   if max((x_vec(loop_x)-GeomPlate(:,1)).^2 + ...
           (y_vec(loop_y)-GeomPlate(:,2)).^2 < (rads.^2).')
    FS_mesh(loop_y,loop_x)=inf;
   end
  end
 end    
end
    
%surf(x_vec,y_vec,FS_mesh)

displ_fs = zeros(y_res,x_res);

th_res = 51; r_res = 31;

th_vec(:,1) = linspace(-pi,pi,th_res); 

r_vec=zeros(Geom_Vec(2),r_res);
for loop_p=1:Geom_Vec(2)
 r_vec(loop_p,:) = linspace(0, rads(loop_p), r_res);
end

displ_ice = zeros(th_res,r_res,Geom_Vec(2));


[Rm,Tm,Rp,Tp,v_vec,u_vec] = ... %displ_fs, displ_ice, 
    fn_MultiMode_MultiFloe(parameter_vector, Vert_Dim, evs, Geom_Vec, ...
    kappa, thicks, rads, [cx;cy], ...
    r_vec, th_vec, x_vec, y_vec, FS_mesh, res, extra_pts);

displ_fs = conj(displ_fs);
displ_ice = conj(displ_ice);

if 0
    figure; hold on
          for loop_P=1:Geom_Vec(2)
           XX = cos(th_vec)*r_vec(loop_P,:); YY = sin(th_vec)*r_vec(loop_P,:);
           surf(cx(loop_P)+XX,cy(loop_P)+YY,real(displ_ice(:,:,loop_P)))
          end
    surf(x_vec,y_vec,real(displ_fs))
    shading interp
    xlabel('x'); ylabel('y'); zlabel('z')
    colorbar

end

% display(['Transmitted amplitudes: ' num2str(abs(Ap(:,1).'))])
% for loop_Y=2:length(Ap(1,:))
%  display(['                        ' num2str(abs(Ap(:,loop_Y).'))])
% end

disp('%----------- END: Steady_MultiFloe ------------%')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Subfns - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Out, D, R, x, y] = fn_GeomDef3d(Geom)

% [Out, D, R, x, y] = GeomDef3d(Geom)
% - Geom=[centre (x,y), radius, thickness] for all plates
% - Out = [scaling No.plates depth thickness l w] (scaled)
% - where l is length and wis width/2 (of wavetank)
% - x and y define the centres of plates

%% Geometry along the x-axis
Np=size(Geom,1);            % Number of plates
x=zeros(1,Np);              % x-coords centre of plates
l=10;                       % length of wavetank
for p=1:Np
    x(p)=Geom(p,1); 
end

%% Geometry along the y-axis
y=zeros(1,Np);              % y-coords centre of plates
w=5;               % walls of wavetank
for p=1:Np
    y(p)=Geom(p,2); 
end

%% Geometry along the z-axis
H=2;                         % in metres

%% Scaling length
L_scale=1; %H;

R = Geom(:,3).'/L_scale; % -scaled radii
D = Geom(:,4).'/L_scale; % -scaled thicknesses
x = x/L_scale; y = y/L_scale; 

%% Output vector of tank parameters 
Out=[L_scale Np H/L_scale l/L_scale, w/L_scale];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
