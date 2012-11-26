function overlap(GeomPlate, TankDim)

% overlap(GeomPlate, TankDim)
% 
% Check that plates don't overlap &/or touch side walls

%% Extract parameters
Np = size(GeomPlate,1);     % Number of disks
cx = GeomPlate(:,1);        % Vector of disk centre x-coord 
cy = GeomPlate(:,2);        % Vector of disk centre y-coord 
rads = GeomPlate(:,3);      % Vector of disk radii


%% Plate/Plate overlap
distij = zeros(Np);
angij = zeros(Np);

for loop1 = 1:Np
    for loop2 = 1:Np
        [distij(loop1,loop2),angij(loop1,loop2)] = coordinate_change(...
            [cx(loop1),cy(loop1)], [cx(loop2),cy(loop2)]);
    end
end

for loop1 = 1:Np-1
    for loop2 = loop1+1:Np
        sumRads = rads(loop1)+rads(loop2);
        if distij(loop1,loop2) < sumRads
            disp(['plates ',num2str(loop1),' & ',...
                num2str(loop2),' overlapping!!!!!'])
            return
        end
        clear sumRads distij
    end
end

%% Plate/Walls overlap
if min(cy-rads) < 0
    [mny,pl] = min(cy-rads);
    disp(['plate ',num2str(pl),' overlaps y=0'])
end

if max(cy-rads) > TankDim(2)
    [mxy,pl] = max(cy-rads);
    disp(['plate ',num2str(pl),' overlaps y=w'])
end

if min(cx-rads) < 0
    [mnx,pl] = min(cx-rads);
    disp(['plate ',num2str(pl),' overlaps x=0'])
end

if max(cx-rads) > TankDim(1)
    [mxx,pl] = max(cx-rads);
    disp(['plate ',num2str(pl),' overlaps x=l'])
end