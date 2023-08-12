function area = get_ele_area(node_ele)
vec1 = node_ele(:, 2) - node_ele(:, 1);
vec2 = node_ele(:, 3) - node_ele(:, 1);

normal = cross(vec1, vec2);

area = 0.5 * norm(normal);

end

