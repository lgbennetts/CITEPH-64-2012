function Mesh = Mesh_FS_def(GeomDisks, TankDim)

% Mesh = Mesh_FS_def(GeomDisks, TankDim)
%
% Definition of the mesh points for the horizontal free-surface region.

%%% Resolution (in points per metres)
res_x = 11;     % along x-axis
res_y = 11;     % along y-axis
res_r = 16;     % along radius of disk

%%% Horizontal space vectors
Mesh.x_vec = linspace(0,TankDim(1),floor(res_x*TankDim(1)));
Mesh.y_vec = linspace(0,TankDim(2),floor(res_y*TankDim(2)));
Mesh.r_vec = cell(1,size(GeomDisks,1));
for loop_p = 1:size(GeomDisks,1)
    Mesh.r_vec{loop_p} = linspace(0, GeomDisks(loop_p,3), ...
        res_r*GeomDisks(loop_p,3));
end
Mesh.th_vec = linspace(-pi,pi,51);

%%% The x vals that bound the plate-covered region
minx0 = min(GeomDisks(:,1) - GeomDisks(:,3));
maxx0 = max(GeomDisks(:,1) + GeomDisks(:,3));

%%% Free-surface grid
Mesh.FS_mesh = zeros(length(Mesh.y_vec),length(Mesh.x_vec));

for loop_x = 1:length(Mesh.x_vec)
    if Mesh.x_vec(loop_x) > minx0 && Mesh.x_vec(loop_x) < maxx0
        for loop_y = 1:length(Mesh.y_vec)
            if max((Mesh.x_vec(loop_x) - GeomDisks(:,1)).^2 + ...
                    (Mesh.y_vec(loop_y) - GeomDisks(:,2)).^2 < ...
                    (GeomDisks(:,3).^2).')
                Mesh.FS_mesh(loop_y, loop_x) = inf;
            end
        end
    end
end