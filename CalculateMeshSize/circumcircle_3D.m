function diameter = circumcircle_3D(triangle)
% To obtain the diameter of the circumcircle of a triangle in 3D Cartesian
% frame.

% shift this triangle and make p1 = [0,0,0].
new_points = Shift_p1toO(triangle);

side1 = new_points(:, 2) - new_points(:, 1);
side2 = new_points(:, 3) - new_points(:, 1);
normal = cross(side1, side2);

% Solve the rotation mapping to plane-xOy.
Rotation_mapping = Rotation_xOy(normal);

% Rotate the triangle with p1 fixed.
new_triangle = Rotation_mapping * new_points;

[r, cn] = circumcircle(new_triangle(1:2, 1:3), 0);
% This function should be download from "MATLAB - APPS - Get More App".
diameter = 2 * r;

return;
end

% EOF