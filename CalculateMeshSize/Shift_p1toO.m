function new_position = Shift_p1toO(Node_ele)
% Shift the triangle in 3d space, let p1 = (0, 0, 0)'.
% (For the rotation that p1 is fixed.)

p1 = Node_ele(:, 1);

new_position = zeros(size(Node_ele, 1), size(Node_ele, 2));

for jj = 1 : size(Node_ele, 2)
    new_position(:, jj) = Node_ele(:, jj) - p1;
end

end

% EOF
