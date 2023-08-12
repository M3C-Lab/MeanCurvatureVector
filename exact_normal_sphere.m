function Normal = exact_normal_sphere(x, y, z)
% To generate the exact solution of Sphere

i_component = 2 * x;
j_component = 2 * y;
k_component = 2 * z;

Normal = [i_component, j_component, k_component]';

% If outside normal is wanted.
Normal = Normal / norm(Normal);

% If mean curvature vector is wanted.
Normal = Normal / 2;
% Norm = 2 / Radius, our radius = 4.

return;
end

% EOF

