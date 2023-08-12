function [normal,Jacobian] = get_normal(Elem_degree, nodes, area_coor)
% To obtain the pointwise normal vector and Jacobian of a quadrature point.
% Input:
%   degree: The degree of triangular element, 1 ~ 2
%   nodes: The nodes of the element.
%   area_coor: The area coordinates of the quadrature point.
% Output:
%   normal: The outside normal at that quadrature point.
%   Jacobian: Point-wise Jacobian at that quadrature point.

a = zeros(3, 1);
b = zeros(3, 1);

for phys_comp = 1 : 3
    a(phys_comp) = Grad(Elem_degree, phys_comp, 1, nodes, area_coor);
    b(phys_comp) = Grad(Elem_degree, phys_comp, 2, nodes, area_coor);
end

normal = cross(a, b);
Jacobian = norm(normal);
% Normalization.
normal = normal / Jacobian;

return;

end

% EOF
