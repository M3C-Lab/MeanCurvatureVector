function new_vector = rotate_xaxis(vector, phi)

rotation = [1, 0, 0;
            0, cos(phi), -sin(phi);
            0, sin(phi), cos(phi)];

new_vector = rotation * vector;
end

