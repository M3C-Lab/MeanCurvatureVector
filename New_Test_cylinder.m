clc;clear all;close all;
% Here is a test for the potential influcence of the number of triangles
% surrounding a node.
addpath("Assembly&Quadrature","CalculateMeshSize", "Preprocess", "ForTest");

msh = load_gmsh2('Cylinder.msh');

Elem_degree = 1;
% The degree of the element.

IEN_v = get_IEN_v(msh, Elem_degree); 
% The IEN of volume-element.

IEN_s = get_IEN_s(msh, Elem_degree, 'surf');
% The IEN of surface-element on the outside surface.

Nodes = unique(IEN_s);
nNode = length(Nodes);

ID = zeros(3, nNode);
for jj = 1 : nNode
    for ii = 1 : 3
        ID(ii, jj) = 3 * (jj - 1) + ii; 
    end
end

[BP_inlet, normal_inlet] = find_BP(msh, Elem_degree, 'inlet');
mcn_inlet = normal_inlet;

[BP_outlet, normal_outlet] = find_BP(msh, Elem_degree, 'outlet');
mcn_outlet = normal_outlet;

node_area = zeros(nNode, 3);
% Here node_area records the surrounding (1-ring) triangles:
% node_area(ii, *): The ii-th node.
% node_area(ii, 1): The total area of surrounding triangles.
% node_area(ii, 2): The total number of surrounding triangles.
% node_area(ii, 3): To show whether the node is on the boundary:
% node_area(ii, 3) = 0 -- the node is on the interior surface;
% node_area(ii, 3) = 1 -- on the 'inlet' boundary;
% node_area(ii, 3) = 2 -- on the 'outlet' boundary;

% Dirichlet BC
Diri_normal = zeros(3, nNode);
Diri_mcn = zeros(3, nNode);
for nn = 1 : length(BP_inlet)
    Diri_normal(:, BP_inlet(nn)) = normal_inlet(:, nn);
    Diri_mcn(:, BP_inlet(nn)) = mcn_inlet(:, nn);
    node_area(BP_inlet(nn), 3) = 1;
end
for nn = 1 : length(BP_outlet)
    Diri_normal(:, BP_outlet(nn)) = normal_outlet(:, nn);
    Diri_mcn(:, BP_outlet(nn)) = mcn_outlet(:, nn);
    node_area(BP_outlet(nn), 3) = 2;
end

% Modify ID array
ID = modify_ID(ID, BP_inlet);
ID = modify_ID(ID, BP_outlet);

LM = get_LM(ID, IEN_s, Nodes);

[tqp, wtqp, ntqp] = TriQuad(6);
[qp, wqp]  = Gauss(6, 0, 1);
% The quadrature rule.

n_eq = 3 * (nNode - length(BP_outlet) - length(BP_inlet));
% The number of equations.
n_EN = size(IEN_s, 1);
% The number of nodes in an element.
n_Elem = size(IEN_s, 2);
% The number of elements.

M = zeros(n_eq, n_eq);

b = zeros(n_eq, 1);

for ee = 1 : n_Elem
    node_ele = zeros(3, n_EN);
    for jj = 1 : n_EN
        node_ele(:, jj) = ...
            [msh.POS(IEN_s(jj, ee), 1), msh.POS(IEN_s(jj, ee), 2), msh.POS(IEN_s(jj, ee), 3)]';
    end


    if Elem_degree == 1
        area_ele = get_ele_area(node_ele);
        for aa = 1 : n_EN
            node_area(IEN_s(aa, ee), 1) = node_area(IEN_s(aa, ee), 1) + area_ele;
            node_area(IEN_s(aa, ee), 2) = node_area(IEN_s(aa, ee), 2) + 1;
        end
    end

    M_ele = zeros(3 * n_EN, 3 * n_EN);

    b_ele = zeros(3 * n_EN, 1);

    for qua = 1 : ntqp
        Phi_matrix = zeros(3, 3 * n_EN);

        for aa = 1 : n_EN           
            Na = TriBasis(Elem_degree, aa, 0, 0, tqp(:, qua));
            Phi_matrix(1 : 3, 3 * aa - 2 : 3 * aa) = [Na,  0,  0;
                                                       0, Na,  0;
                                                       0,  0, Na];
        end

        % (Normal and Jacobian)
        [normal_f, Jacobian] = get_normal(Elem_degree, node_ele, tqp(:, qua));

        M_ele = M_ele + Jacobian * wtqp(qua) * (Phi_matrix' * Phi_matrix);

        b_ele = b_ele + Jacobian * wtqp(qua) * (Phi_matrix' * -normal_f);
    end

    for aa = 1 : 3 * n_EN
        LM_a = LM(aa, ee);
        if LM_a > 0
            b(LM_a) = b(LM_a) + b_ele(aa);
            for bb = 1 : 3 * n_EN
                LM_b = LM(bb, ee);
                if LM_b > 0
                    M(LM_a, LM_b) = M(LM_a, LM_b) + M_ele(aa, bb);
                else
                    % Dirichlet BC
                    bb_node = floor((bb-1) / 3) + 1;
                    g_normal = Diri_normal(:, IEN_s(bb_node, ee));
                    bb_direc = mod(bb, 3);
                    if bb_direc == 0
                        bb_direc = 3;
                    end
                    b(LM_a) = b(LM_a) - M_ele(aa, bb) * g_normal(bb_direc);
                end
            end
        end
    end
end

nh_normal = M \ b;

result_normal = zeros(3, nNode);

for ii = 1 : 3
    for jj = 1 : nNode
        if ID(ii, jj) ~= 0
            result_normal(ii, jj) = nh_normal(ID(ii, jj));
        end
    end
end

result_normal = result_normal + Diri_normal;

node_proj = cell(nNode, 1);
for nn = 1 : nNode
    I = eye(3);
    tensor_nn = result_normal(:, nn) * result_normal(:, nn)';
    node_proj{nn, 1} = I - tensor_nn;
end


M = zeros(n_eq, n_eq);

Sx = zeros(n_eq, 1);

hh = 0.0;

for ee = 1 : n_Elem
    node_ele = zeros(3, n_EN);
    for jj = 1 : n_EN
        node_ele(:, jj) = ...
            [msh.POS(IEN_s(jj, ee), 1), msh.POS(IEN_s(jj, ee), 2), msh.POS(IEN_s(jj, ee), 3)]';
    end

    triangle = node_ele(:, 1:3);
    hh_ele = circumcircle_3D(triangle);
    if hh_ele > hh
        hh = hh_ele;
    end

    M_ele = zeros(3 * n_EN, 3 * n_EN);

    Sx_ele = zeros(3 * n_EN, 1);

    for qua = 1 : ntqp
        Phi_matrix = zeros(3, 3 * n_EN);

        Sx_qua = zeros(3 * n_EN, 1);

        if Elem_degree == 1
            dx_proj = zeros(3, 3);
            % Approximate dx / dx Projection

            for aa = 1 : n_EN
                dx_proj = dx_proj + node_proj{IEN_s(aa, ee), 1} * ...
                    TriBasis(Elem_degree, aa, 0, 0, tqp(:,qua));
            end
        end
        
        F_hat = zeros(3, 2);
        for ii = 1 : 3
            F_hat(ii, :) = [Grad(Elem_degree, ii, 1, node_ele, tqp(:,qua)),...
                            Grad(Elem_degree, ii, 2, node_ele, tqp(:,qua))];
        end

        G_FFF = F_hat' * F_hat;

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

                    if Elem_degree ~= 1
                        vec2(nn, 1) = Xn_vec / G_FFF * Xm_vec;
                    else
                        vec2(nn, 1) = dx_proj(mm, nn);
                    end
                end
                Sx_qua(3 * aa - 3 + mm, 1) = dot(vec1, vec2);
            end
        end

        % (Normal and Jacobian)
        [normal_f, Jacobian] = get_normal(Elem_degree, node_ele, tqp(:, qua));

        M_ele = M_ele + Jacobian * wtqp(qua) * (Phi_matrix' * Phi_matrix);

        Sx_ele = Sx_ele + Jacobian * wtqp(qua) * Sx_qua;
    end

    for aa = 1 : 3 * n_EN
        LM_a = LM(aa, ee);
        if LM_a > 0
            Sx(LM_a) = Sx(LM_a) + Sx_ele(aa);
            for bb = 1 : 3 * n_EN
                LM_b = LM(bb, ee);
                if LM_b > 0
                    M(LM_a, LM_b) = M(LM_a, LM_b) + M_ele(aa, bb);
                else
                    % Dirichlet BC
                    bb_node = floor((bb-1) / 3) + 1;
                    g_mcn = Diri_mcn(:, IEN_s(bb_node, ee));
                    bb_direc = mod(bb, 3);
                    if bb_direc == 0
                        bb_direc = 3;
                    end
                    Sx(LM_a) = Sx(LM_a) - M_ele(aa, bb) * g_mcn(bb_direc);
                end
            end
        end
    end
end

nh_mcn = M \ Sx;

result_mcn = zeros(3, nNode);

for ii = 1 : 3
    for jj = 1 : nNode
        if ID(ii, jj) ~= 0
            result_mcn(ii, jj) = nh_mcn(ID(ii, jj));
        end
    end
end

result_mcn = result_mcn + Diri_mcn;


result = result_mcn;

X = zeros(nNode, 1);
Y = zeros(nNode, 1);
Z = zeros(nNode, 1);
for jj = 1 : nNode
    X(jj) = msh.POS(Nodes(jj), 1);
    Y(jj) = msh.POS(Nodes(jj), 2);
    Z(jj) = msh.POS(Nodes(jj), 3);
end
U = result(1, :)';
V = result(2, :)';
W = result(3, :)';

MESH = alphaShape(X, Y, Z, 4, 'HoleThreshold', 1e-6);

figure(1);
quiver3(X, Y, Z, U, V, W);
hold on;
plot(MESH);

for ii = 1 : nNode
    vector = result(:, ii);
    fprintf('number:%d,  norm: %f\n\n', ii, norm(vector));
end

% Record the minimum and maximum norm:
max_norm = 0.0;
min_norm = 10000.0;
for ii = 1 : nNode
    vector = result(:, ii);
    % fprintf('number:%d,  norm: %f\n\n', ii, norm(vector));
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