function [boundary_points, normals_BC] = find_BP(msh, Elem_degree, boundary_name)
% To find the points on the boundary of the open surface, and generate
% their normals

for I = 1 : size(msh.PhyGrp, 1)
    if strcmp(boundary_name, msh.PhyGrp{I, 3})
        PG_number = msh.PhyGrp{I, 2};
        break;
    end
end

boundary_points = [ ];

if Elem_degree == 1
    for ii = 1 : msh.nbLines
        if msh.LINES(ii, 3) == PG_number
            nodes = msh.LINES(ii, 1:2);
            boundary_points = [boundary_points, nodes];
        end
    end
elseif Elem_degree == 2
    for ii = 1 : msh.nbLines3
        if msh.LINES3(ii, 4) == PG_number
            nodes = msh.LINES3(ii, 1:3);
            boundary_points = [boundary_points, nodes];
        end
    end
end

boundary_points = unique(boundary_points);

total_x = 0;
total_y = 0;
total_z = 0;

for jj = 1 : length(boundary_points)
    total_x = total_x + msh.POS(boundary_points(jj), 1);
    total_y = total_y + msh.POS(boundary_points(jj), 2);
    total_z = total_z + msh.POS(boundary_points(jj), 3);
end

average_x = total_x / length(boundary_points);
average_y = total_y / length(boundary_points);
average_z = total_z / length(boundary_points);

normals_BC = zeros(3, length(boundary_points));

for kk = 1 : length(boundary_points)
    normal = [msh.POS(boundary_points(kk), 1) - average_x;
              msh.POS(boundary_points(kk), 2) - average_y;
              msh.POS(boundary_points(kk), 3) - average_z];
    normal = normal / norm(normal);
    normals_BC(1:3, kk) = normal;
end

return;

end

% EOF
