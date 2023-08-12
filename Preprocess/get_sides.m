function [triangle_side, side_node] = get_sides(Elem_degree, IEN, nElem_l)
% To get the triangle-side information and the side-node information.
% Input:
%   nElem_l: The total number of line elements on the boundary of open surface.
%   IEN: Surface IEN.
% Output:
%   side_node: The side-node relationship.
%       side_node(*, ii): the side ii, 
%       side_node(nn, ii):the node nn of side ii.

%   triangle_side: The triangle-side relationship.
%       triangle_side(*, ee): the triangle element ee,
%       triangle_side(nn, ee): the side nn of element ee.
%
% positive index means the side is identified counterclockwise, 
% and negative index means clockwise.
% Example:
%           node3  -------- node 4
%                 /\      /
%                /  \    /
% tri-element 1 /    \  /  tri-element 2
%              /      \/
%              -------
%          node 1     node2
%
%  (The outside normal points to the outside of the screen.)
%
%  If we record "side 56:[node 2; node 3]" in side_node, i.e. side_node(:, 56) = [2;3],
%  then it will be noted as "56" for element 1, but "-56" for element 2 in triangle_side.


    n_Elem = size(IEN, 2);
    max_n_side = 1.5 * n_Elem + 0.5 * nElem_l;

    triangle_side = zeros(3, n_Elem);

    if Elem_degree == 1
        side_node = zeros(2, max_n_side);
    end

    side_number = 1;

    for ee = 1 : n_Elem
        if Elem_degree == 1
            sides = [IEN(1, ee), IEN(2, ee); IEN(2, ee), IEN(3, ee); IEN(3, ee), IEN(1, ee)]';
        end
        
        for ss = 1 : 3 
            for nn = 1 : side_number
                if side_node(1, nn) == sides(1, ss) && side_node(2, nn) == sides(2, ss)
                    triangle_side(ss, ee) = nn;
                    break;
                elseif side_node(1, nn) == sides(2, ss) && side_node(2, nn) == sides(1, ss)
                    triangle_side(ss, ee) = -nn;
                    break;
                end

                if nn == side_number    % If there is no recorded side information.
                    side_node(:, side_number) = sides(:, ss);
                    triangle_side(ss, ee) = nn;
                    side_number  = side_number + 1;     % Add the new side.
                end
            end
        end
    end

    return;

end

% EOF

