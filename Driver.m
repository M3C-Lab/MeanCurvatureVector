clc;clear all;close all;
% function 'circumcircle' should be download from "MATLAB - APPS - Get More App".
addpath("Assembly&Quadrature","CalculateMeshSize", "Preprocess");

msh = load_gmsh2('Sphere.msh');

Elem_degree = 1;
% The degree of the element.

gamma = 0.0001;

IEN_v = get_IEN_v(msh, Elem_degree); 
% The IEN of volume-element.

IEN_s = get_IEN_s(msh, Elem_degree, 'surf');
% The IEN of surface-element on the outside surface.

Nodes = unique(IEN_s);
% The nodes

nNode = length(Nodes);  
% The number of node on the outside surface.

ID = zeros(3, nNode);
for jj = 1 : nNode
    for ii = 1 : 3
        ID(ii, jj) = 3 * (jj - 1) + ii; 
    end
end
LM = get_LM(ID, IEN_s, Nodes);


if Elem_degree == 1
    facet_normal = get_facet_normal(msh, IEN_v, IEN_s); 
    % The facet normal on the outside surface.

    [tri_side, side_node] = get_sides(Elem_degree, IEN_s, msh.nbLines);
    % The triangle-side and side-node relationship.

    find_node_equa = zeros(6, size(side_node, 2));
    for ii = 1 : 2
        for jj = 1 : size(side_node, 2)
            for kk = 1 : 3
                find_node_equa(3 * ii - 3 + kk, jj) = 3 * side_node(ii, jj) - 3 + kk;
            end
        end
    end
    % The 'LM' array for side_node array.

    side_tri = side_triangle(size(side_node, 2), tri_side);
    % The side-triangle relationship.

    conormal = get_conormal(facet_normal, side_tri, side_node, msh);
    % The side-conormal relationship
end

[tqp, wtqp, ntqp] = TriQuad(6);
% The quadrature rule for triangle.

[qp, wqp]  = Gauss(6, 0, 1);
% The quadrature rule for line element.

n_eq = 3 * nNode;
% The number of equations.
n_EN = size(IEN_s, 1);
% The number of nodes in an element.
n_Elem = size(IEN_s, 2);
% The number of elements.

M = sparse(n_eq, n_eq);

b = zeros(n_eq, 1);

Sx = zeros(n_eq, 1);

% For the edge stablization we have to calculate mesh size in advance.
hh = 0.0;

%Assembly.
tic;

% Record the 1-ring element number of the node (for 1st order element).
node_element_number = zeros(nNode, 1);

for ee = 1 : n_Elem
    node_ele = zeros(3, n_EN);
    for jj = 1 : n_EN
        node_ele(:, jj) = ...
            [msh.POS(IEN_s(jj, ee), 1), msh.POS(IEN_s(jj, ee), 2), msh.POS(IEN_s(jj, ee), 3)]';
        node_element_number(IEN_s(jj, ee), 1) = node_element_number(IEN_s(jj, ee), 1) + 1;
    end
    
    % For the edge stablization we have to calculate mesh size in advance.
    hh = 0.0;
    triangle = node_ele(:, 1:3);
    hh_ele = circumcircle_3D(triangle);
    if hh_ele > hh
        hh = hh_ele;
    end
    % The mesh size of paper is the diameter of the largest inscribed
    % circle of elements, but here we use circumcircle. 

    M_ele = zeros(3 * n_EN, 3 * n_EN);

    b_ele = zeros(3 * n_EN, 1);

    Sx_ele = zeros(3 * n_EN, 1);

    for qua = 1 : ntqp
        Phi_matrix = zeros(3, 3 * n_EN);

        Sx_qua = zeros(3 * n_EN, 1);
        
        F_hat = zeros(3, 2);
        for ii = 1 : 3
            F_hat(ii, :) = [Grad(Elem_degree, ii, 1, node_ele, tqp(:,qua)),...
                            Grad(Elem_degree, ii, 2, node_ele, tqp(:,qua))];
        end
        % The 'deformation gradient' of isoparametric mapping.

        G_FFF = F_hat' * F_hat;
        % The First Fundamental Form.

        for aa = 1 : n_EN           
            Na = TriBasis(Elem_degree, aa, 0, 0, tqp(:, qua));
            Phi_matrix(1 : 3, 3 * aa - 2 : 3 * aa) = [Na,  0,  0;
                                                       0, Na,  0;
                                                       0,  0, Na];

            % (For mean curvature vector)
            Na_xi_vec = [TriBasis(Elem_degree, aa, 1, 0, tqp(:,qua));
                     TriBasis(Elem_degree, aa, 0, 1, tqp(:,qua))];
            for mm = 1 : 3
                vec1 = zeros(3, 1);
                vec2 = zeros(3, 1);
                Xm_vec = F_hat(mm, :)';
                for nn = 1 : 3
                    Xn_vec = F_hat(nn, :);
                    vec1(nn, 1) = Xn_vec / G_FFF * Na_xi_vec;
                    vec2(nn, 1) = Xn_vec / G_FFF * Xm_vec;
                end
                Sx_qua(3 * aa - 3 + mm, 1) = dot(vec1, vec2);
            end
        end

        % (Normal and Jacobian)
        [normal_f, Jacobian] = get_normal(Elem_degree, node_ele, tqp(:, qua));

        M_ele = M_ele + Jacobian * wtqp(qua) * (Phi_matrix' * Phi_matrix);

        b_ele = b_ele + Jacobian * wtqp(qua) * (Phi_matrix' * normal_f);

        Sx_ele = Sx_ele + Jacobian * wtqp(qua) * Sx_qua;
    end

    for aa = 1 : 3 * n_EN
        LM_a = LM(aa, ee);
        if LM_a > 0
            b(LM_a) = b(LM_a) + b_ele(aa);
            Sx(LM_a) = Sx(LM_a) + Sx_ele(aa);
            for bb = 1 : 3 * n_EN
                LM_b = LM(bb, ee);
                if LM_b > 0
                    M(LM_a, LM_b) = M(LM_a, LM_b) + M_ele(aa, bb);
                end
            end
        end
    end
end
toc;

tic;

% Edge stablization
if Elem_degree == 1
    for ss = 1 : size(side_tri, 2)
        if  side_tri(2, ss) ~= 0 % Not on the boundary.
            
            % Left element:
            ele_K1 = abs(side_tri(1, ss));
            [node_number_K1, new_qp_K1] = find_sidenode_number(side_node(:, ss), IEN_s(:, ele_K1), qp);
            conormal_K1 = conormal(1:3, ss);

            % Right element:
            ele_K2 = abs(side_tri(2, ss));
            [node_number_K2, new_qp_K2] = find_sidenode_number(side_node(:, ss), IEN_s(:, ele_K2), qp);
            conormal_K2 = conormal(4:6, ss);

            node_ele_K1 = zeros(3, n_EN);
            node_ele_K2 = zeros(3, n_EN);
            for jj = 1 : n_EN
                node_ele_K1(:, jj) = ...
                    [msh.POS(IEN_s(jj, ele_K1), 1), msh.POS(IEN_s(jj, ele_K1), 2), msh.POS(IEN_s(jj, ele_K1), 3)]';

                node_ele_K2(:, jj) = ...
                    [msh.POS(IEN_s(jj, ele_K2), 1), msh.POS(IEN_s(jj, ele_K2), 2), msh.POS(IEN_s(jj, ele_K2), 3)]';
            end

            side_vector = [msh.POS(side_node(2, ss), 1) - msh.POS(side_node(1, ss), 1);
                           msh.POS(side_node(2, ss), 2) - msh.POS(side_node(1, ss), 2);
                           msh.POS(side_node(2, ss), 3) - msh.POS(side_node(1, ss), 3)];
            
            dx_dxi = norm(side_vector); % The Jacobian for line quadrature.

            J_side = zeros(2, 2);

            for qua = 1 : length(qp)

                F_hat_K1 = zeros(3, 2);
                F_hat_K2 = zeros(3, 2);
                for ii = 1 : 3
                    F_hat_K1(ii, :) = [Grad(Elem_degree, ii, 1, node_ele_K1, new_qp_K1(:,qua)),...
                                       Grad(Elem_degree, ii, 2, node_ele_K1, new_qp_K1(:,qua))];

                    F_hat_K2(ii, :) = [Grad(Elem_degree, ii, 1, node_ele_K2, new_qp_K2(:,qua)),...
                                       Grad(Elem_degree, ii, 2, node_ele_K2, new_qp_K2(:,qua))];
                end

                G_FFF_K1 = F_hat_K1' * F_hat_K1;
                G_FFF_K2 = F_hat_K2' * F_hat_K2;

                F_co_K1 = F_hat_K1' * conormal_K1;
                F_co_K2 = F_hat_K2' * conormal_K2;

                N1_grad_K1 = [TriBasis(Elem_degree, node_number_K1(1), 1, 0, new_qp_K1(:,qua));
                              TriBasis(Elem_degree, node_number_K1(1), 0, 1, new_qp_K1(:,qua))];

                N1_grad_K2 = [TriBasis(Elem_degree, node_number_K2(1), 1, 0, new_qp_K2(:,qua));
                              TriBasis(Elem_degree, node_number_K2(1), 0, 1, new_qp_K2(:,qua))];

                N2_grad_K1 = [TriBasis(Elem_degree, node_number_K1(2), 1, 0, new_qp_K1(:,qua));
                              TriBasis(Elem_degree, node_number_K1(2), 0, 1, new_qp_K1(:,qua))];

                N2_grad_K2 = [TriBasis(Elem_degree, node_number_K2(2), 1, 0, new_qp_K2(:,qua));
                              TriBasis(Elem_degree, node_number_K2(2), 0, 1, new_qp_K2(:,qua))];

                N1K1 = dot(F_co_K1', inv(G_FFF_K1)* N1_grad_K1); 
                N1K2 = dot(F_co_K2', inv(G_FFF_K2) * N1_grad_K2);
                N2K1 = dot(F_co_K1', inv(G_FFF_K1) * N2_grad_K1); 
                N2K2 = dot(F_co_K2', inv(G_FFF_K2) * N2_grad_K2);

                N1_term = (N1K1 + N1K2);
                N2_term = (N2K1 + N2K2);

                N1N1 = N1_term^2;
                N1N2 = N1_term * N2_term;
                N2N2 = N2_term^2;

                J_qua = [N1N1, N1N2;
                         N1N2, N2N2];

                J_side = J_side + J_qua * wqp(qua) * dx_dxi;

            end
            
            % Additional assembly
            for aa = 1 : 3
                M(find_node_equa(aa, ss), find_node_equa(aa, ss)) = ...
                    M(find_node_equa(aa, ss), find_node_equa(aa, ss)) + gamma * hh * J_side(1, 1);

                M(find_node_equa(aa, ss), find_node_equa(aa + 3, ss)) = ...
                    M(find_node_equa(aa, ss), find_node_equa(aa + 3, ss)) + gamma * hh * J_side(1, 2);

                M(find_node_equa(aa + 3, ss), find_node_equa(aa, ss)) = ...
                    M(find_node_equa(aa + 3, ss), find_node_equa(aa, ss)) + gamma * hh * J_side(2, 1);

                M(find_node_equa(aa + 3, ss), find_node_equa(aa + 3, ss)) = ...
                    M(find_node_equa(aa + 3, ss), find_node_equa(aa + 3, ss)) + gamma * hh * J_side(2, 2);
            end
        end
    end
end
toc;

% % % To get outside normal.
% % nh = M \ b;

% To get mean curvature vector.
nh = M \ Sx;


nh_result = zeros(3, nNode);
for ii = 1 : 3
    for jj = 1 : nNode
        nh_result(ii, jj) = nh(ID(ii, jj));
    end
end
% nh_result(:, ii): The vector result of the ii-th node.

% Plot
X = zeros(nNode, 1);
Y = zeros(nNode, 1);
Z = zeros(nNode, 1);
NORM = zeros(nNode, 1);
for jj = 1 : nNode
    X(jj) = msh.POS(Nodes(jj), 1);
    Y(jj) = msh.POS(Nodes(jj), 2);
    Z(jj) = msh.POS(Nodes(jj), 3);
    NORM(jj) = norm(nh_result(:, jj));
end
U = nh_result(1, :)';
V = nh_result(2, :)';
W = nh_result(3, :)';

MESH = alphaShape(X, Y, Z, 4, 'HoleThreshold', 1e-6);

figure(1);

quiver3(X, Y, Z, U, V, W);

hold on;
plot(MESH);

figure(2);
scatter3(X,Y,Z, 100, NORM, 'filled');
colorbar;

% Print the norm of each vector.
% Record the minimum and maximum norm.
max_norm = 0.0;
min_norm = 10000.0;
for ii = 1 : nNode
    vector = nh_result(:, ii);
    fprintf('number:%d,  norm: %f\n\n', ii, norm(vector));
    if norm(vector) > max_norm
        max_norm = norm(vector);
        max_norm_node = ii;
    end
    if norm(vector) < min_norm
        min_norm  = norm(vector);
        min_norm_node = ii;
    end
end

fprintf('Min_norm: %f Node: %d;  Max_norm: %f Node: %d\n\n', min_norm, min_norm_node, max_norm, max_norm_node);


% Calculate the relative error:
[tqp, wtqp, ntqp] = TriQuad(6);

top = 0.0; bot = 0.0;
hh = 0.0; % Mesh size.

for ee = 1 : n_Elem
    node_ele = zeros(3, n_EN);
    nh_ele = zeros(3, n_EN);
    for aa = 1 : n_EN
        node_ele(:, aa) = ...
            [msh.POS(IEN_s(aa, ee), 1), ...
             msh.POS(IEN_s(aa, ee), 2), ...
             msh.POS(IEN_s(aa, ee), 3)]';
        found_node_number = find(Nodes == IEN_s(aa, ee));
        nh_ele(:, aa) = nh_result(:, found_node_number);
    end

    triangle = node_ele(:, 1:3);
    hh_ele = circumcircle_3D(triangle);
    if hh_ele > hh
        hh = hh_ele;
    end

    for qua = 1 : ntqp
        x_qua = 0.0; y_qua = 0.0; z_qua = 0.0;
        nh_qua = zeros(3, 1);
        for aa = 1 :n_EN
            Na = TriBasis(Elem_degree, aa, 0, 0, tqp(:, qua));
            x_qua = x_qua + node_ele(1, aa) * Na;
            y_qua = y_qua + node_ele(2, aa) * Na;
            z_qua = z_qua + node_ele(3, aa) * Na;
            nh_qua = nh_qua + nh_ele(:, aa) * Na;
        end

        [normal_f, Jacobian] = get_normal(Elem_degree, node_ele, tqp(:, qua));
        
        n_exact = exact_normal_sphere(x_qua, y_qua, z_qua);

        top = top + wtqp(qua)* Jacobian * dot((nh_qua - n_exact), (nh_qua - n_exact));
        bot = bot + wtqp(qua)* Jacobian * dot(n_exact, n_exact);
    end
end

top = sqrt(top); bot = sqrt(bot);
rel_error = top / bot;

fprintf('relative error: %f,  mesh size: %f\n\n, gamma: %f\n\n', rel_error, hh, gamma);

% I think the stablization strategy failed.

% EOF