function rotation_mapping = Rotation_xOy(standard)
% Rotate a triangle in 3d space to the xOy plane,
%   by rotate the standard normal to z+ direction.
% Output:
%   rotation_mapping: the rotation matrix that
%       rotation_mapping * standard_normal = (0, 0, |standard_normal|)'.

ex = [1, 0, 0]';
ey = [0, 1, 0]';
ez = [0, 0, 1]';

if standard == ez * norm(standard)
    rotation_mapping = eye(3);
    return;
elseif standard ==  - ez * norm(standard)
    rotation_mapping = - eye(3);
    return;
else
    sub_normal = cross(ez, standard);
    cos_phi = dot(ex, sub_normal) / norm(sub_normal);
    sin_phi = dot(ey, sub_normal) / norm(sub_normal);
    % theta is the angle between the sub nomral and the x-axis.
    % The rotation axis is the z-axis.
    rotation_z = [cos_phi, sin_phi, 0;
                  -sin_phi, cos_phi, 0;
                  0, 0, 1];
    % The inverse rotation over z-axis using (-phi).
    
    cos_theta = dot(ez, standard) / norm(standard);
    sin_theta = norm(cross(ez, standard)) / norm(standard);
    % theta is the angle between the standard nomral and the z-axis.
    % The rotation axis is the x-axis.
    rotation_x = [1, 0, 0;
                  0, cos_theta, sin_theta;
                  0, -sin_theta, cos_theta];
    % The inverse rotation over x-axis using (-theta).
    
    rotation_mapping = rotation_x * rotation_z;
    return;
end

end

% EOF


