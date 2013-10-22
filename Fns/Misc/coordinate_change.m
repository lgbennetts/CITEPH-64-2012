% [R,psi]=coordinate_change(pos_j,pos_l)
%
% This function is used to express the the location of the center of floe l
% in the local coordinate system of floe j.

function [R,psi]=coordinate_change(pos_j,pos_l)

R=sqrt((pos_j(1)-pos_l(1))^2+(pos_j(2)-pos_l(2))^2);
psi=atan2(pos_l(2)-pos_j(2),pos_l(1)-pos_j(1));

if psi<0
    psi=psi+2*pi;
end

return