function Grad_comp = Grad(Elem_degree, phys_comp, para_comp, nodes, area_coor)
% To get one of the (partial derivative) gradient components of parametric mapping.
% For example: dx/dr.
% Input:
%   degree: The degree of triangular element, 1 ~ 2
%   n_EN: The number of nodes in an element.
%   phys_comp: x -- 1, y -- 2, z -- 3.
%   para-comp: r -- 1, s -- 2.
%   nodes: The nodes of the element.
%   area_coor: The area coordinates of the point.
% Output: 
%   Grad_comp: the value of gradient components.

Grad_comp = 0;

if para_comp == 1
    for ii = 1 : size(nodes, 2)
        Grad_comp = Grad_comp + nodes(phys_comp, ii) * TriBasis(Elem_degree, ii, 1, 0, area_coor);
    end
elseif para_comp == 2
    for ii = 1 : size(nodes, 2)
        Grad_comp = Grad_comp + nodes(phys_comp, ii) * TriBasis(Elem_degree, ii, 0, 1, area_coor);
    end
end

return;

end

% EOF
