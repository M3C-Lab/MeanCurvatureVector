function s_t = side_triangle(n_side, t_s)
% Input:
%   n_side: The total number of the sides.
%   t_s: The triangle-side relationship.
% Output:
%   s_t: The side-triangle relationship.
%       s_t(*, ii): The side ii.
%       s_t(1, ii): The first element of side ii.
%       s_t(2, ii): The second element of side ii.
%           If the side is on the boundary, it will be remained as zero.
%
% positive index means the element is on the left of the side, and
% negative index means the right.
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
%  If we record "side 56:[node 2; node 3]" in side_node,then the s_t will
%  be "side 56:[1, -2]" or "side 56:[-2, 1]"

s_t = zeros(2, n_side);

for ee = 1 : size(t_s, 2)
    for ss = 1 : 3
        if t_s(ss, ee) > 0
            ele = ee;
        elseif t_s(ss, ee) < 0
            ele = -ee;
        end

        if s_t(1, abs(t_s(ss, ee))) == 0
            s_t(1, abs(t_s(ss, ee))) = ele;
        else
            s_t(2, abs(t_s(ss, ee))) = ele;
        end
    end
end

end

% EOF
