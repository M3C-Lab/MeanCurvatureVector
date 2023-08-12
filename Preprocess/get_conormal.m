function conormal = get_conormal(facet_normal, side_tri, side_node, msh)
% To get the conormal list corresponding to 'side_tri' array.
% Input:
%   facet_normal: Facet normal of all triangles.
%   side_tri: The side-triangle relationship.
%   side_node: The side-node relationship.
%   msh: info from msh file.
% Output:
%   conormal:
%       conormal(*, ii): the conormals of the ii-th side:
%       conormal(1:3, ii): the conormal of first triangle of the ii-th side.
%       conormal(4:6, ii): the conormal of second triangle of the ii-th side.


n_side = size(side_tri, 2);
conormal = zeros(6, n_side);

for ii = 1 : n_side
    if side_tri(2, ii) ~= 0
        side_vector = [msh.POS(side_node(2, ii), 1) - msh.POS(side_node(1, ii), 1);
                       msh.POS(side_node(2, ii), 2) - msh.POS(side_node(1, ii), 2);
                       msh.POS(side_node(2, ii), 3) - msh.POS(side_node(1, ii), 3)];

        side_vector = side_vector / norm(side_vector);

        if side_tri(1, ii) > 0
            tK1 = cross(facet_normal(:, side_tri(1, ii)), side_vector);
            tK2 = cross(side_vector, facet_normal(:, -side_tri(2, ii)));
        elseif side_tri(1, ii) < 0
            tK1 = cross(side_vector, facet_normal(:, -side_tri(1, ii)));
            tK2 = cross(facet_normal(:, side_tri(2, ii)), side_vector);
        end


        conormal(1:3, ii) = tK1;
        conormal(4:6, ii) = tK2;
    end
end

end

% EOF

