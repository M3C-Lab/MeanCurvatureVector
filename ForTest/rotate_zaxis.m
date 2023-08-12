function new_vector = rotate_zaxis(vector, theta)

rotation = [cos(theta), -sin(theta), 0;
            sin(theta), cos(theta), 0;
            0, 0, 1];

new_vector = rotation * vector;
end

