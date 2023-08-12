clc; clear all;close all;
% Here is a test for the potential influcence of the number of triangles
% surrounding a node.

addpath("Assembly&Quadrature","CalculateMeshSize", "Preprocess", 'ForTest');

Elem_degree = 1;

gamma1 = 1;

n_EN = 3;

[tqp, wtqp, ntqp] = TriQuad(6);
[qp, wqp]  = Gauss(6, 0, 1);
% The quadrature rule.

p1 = [0, 0, 4]';
% Suppose a point p1 on the top of a sphere of radius = 4

%% Change this
n_Elem = 6;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 element
if n_Elem == 6
    p2 = rotate_xaxis(p1, pi / 15);
    p3 = rotate_zaxis(p2, pi / 3);
    p4 = rotate_zaxis(p3, pi / 3);
    p5 = rotate_zaxis(p4, pi / 3);
    p6 = rotate_zaxis(p5, pi / 3);
    p7 = rotate_zaxis(p6, pi / 3);
    % Construct 6 surrounding elements
    
    nNode = 7;
    Nodes = [1, 2, 3, 4, 5, 6, 7];
    
    POS = [p1, p2, p3, p4, p5, p6, p7];
    
    IEN_s = [1, 1, 1, 1, 1, 1;
             2, 3, 4, 5, 6, 7;
             3, 4, 5, 6, 7, 2];
    
    ID = [1, 0, 0, 0, 0, 0, 0;
          2, 0, 0, 0, 0, 0, 0;
          3, 0, 0, 0, 0, 0, 0];
    
    LM = get_LM(ID, IEN_s, Nodes);
    
    % Design Dirichlet BC for p2~p7.
    Diri_normal = zeros(3, 7);
    Diri_mcn = zeros(3, 7);
    for nn = 2 : 7
        point = POS(:, nn);
        normal = point / norm(point); % Outside normal.
        Diri_normal(:, nn) = normal;

        mcn = normal / 2; % Mean curvature normal = 2 * Outside normal / Radius.
        Diri_mcn(:, nn) = mcn;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 element
if n_Elem == 5
    p2 = rotate_xaxis(p1, pi / 15);
    p3 = rotate_zaxis(p2, 2*pi / 5);
    p4 = rotate_zaxis(p3, 2*pi / 5);
    p5 = rotate_zaxis(p4, 2*pi / 5);
    p6 = rotate_zaxis(p5, 2*pi / 5);
    
    nNode = 6;
    Nodes = [1, 2, 3, 4, 5, 6];
    
    POS = [p1, p2, p3, p4, p5, p6];
    
    IEN_s = [1, 1, 1, 1, 1;
             2, 3, 4, 5, 6;
             3, 4, 5, 6, 2];
    
    ID = [1, 0, 0, 0, 0, 0;
          2, 0, 0, 0, 0, 0;
          3, 0, 0, 0, 0, 0];
    
    LM = get_LM(ID, IEN_s, Nodes);
    
    Diri_normal = zeros(3, 6);
    Diri_mcn = zeros(3, 6);
    for nn = 2 : 6
        point = POS(:, nn);
        normal = point / norm(point);
        Diri_normal(:, nn) = normal;
        mcn = normal / 2;
        Diri_mcn(:, nn) = mcn;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n_Elem == 4
    p2 = rotate_xaxis(p1, pi / 15);
    p3 = rotate_zaxis(p2, pi / 2);
    p4 = rotate_zaxis(p3, pi / 2);
    p5 = rotate_zaxis(p4, pi / 2);
    
    nNode = 5;
    Nodes = [1, 2, 3, 4, 5];
    
    POS = [p1, p2, p3, p4, p5];
    
    IEN_s = [1, 1, 1, 1;
             2, 3, 4, 5;
             3, 4, 5, 2];
    
    ID = [1, 0, 0, 0, 0;
          2, 0, 0, 0, 0;
          3, 0, 0, 0, 0];
    
    LM = get_LM(ID, IEN_s, Nodes);
    
    Diri_normal = zeros(3, 5);
    Diri_mcn = zeros(3, 5);
    for nn = 2 : 5
        point = POS(:, nn);
        normal = point / norm(point);
        Diri_normal(:, nn) = normal;
        mcn = normal / 2;
        Diri_mcn(:, nn) = mcn;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7 element
if n_Elem == 7
    p2 = rotate_xaxis(p1, pi / 15);
    p3 = rotate_zaxis(p2, 2 * pi / 7);
    p4 = rotate_zaxis(p3, 2 * pi / 7);
    p5 = rotate_zaxis(p4, 2 * pi / 7);
    p6 = rotate_zaxis(p5, 2 * pi / 7);
    p7 = rotate_zaxis(p6, 2 * pi / 7);
    p8 = rotate_zaxis(p7, 2 * pi / 7);
    
    nNode = 8;
    Nodes = [1, 2, 3, 4, 5, 6, 7, 8];
    
    POS = [p1, p2, p3, p4, p5, p6, p7, p8];
    
    IEN_s = [1, 1, 1, 1, 1, 1, 1;
             2, 3, 4, 5, 6, 7, 8;
             3, 4, 5, 6, 7, 8, 2];
    
    ID = [1, 0, 0, 0, 0, 0, 0, 0;
          2, 0, 0, 0, 0, 0, 0, 0;
          3, 0, 0, 0, 0, 0, 0, 0];
    
    LM = get_LM(ID, IEN_s, Nodes);
    
    Diri_normal = zeros(3, 8);
    Diri_mcn = zeros(3, 8);
    for nn = 2 : 8
        point = POS(:, nn);
        normal = point / norm(point);
        Diri_normal(:, nn) = normal;
        mcn = normal / 2;
        Diri_mcn(:, nn) = mcn;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 elements
if n_Elem == 8
    p2 = rotate_xaxis(p1, pi / 15);
    p3 = rotate_zaxis(p2, pi / 4);
    p4 = rotate_zaxis(p3, pi / 4);
    p5 = rotate_zaxis(p4, pi / 4);
    p6 = rotate_zaxis(p5, pi / 4);
    p7 = rotate_zaxis(p6, pi / 4);
    p8 = rotate_zaxis(p7, pi / 4);
    p9 = rotate_zaxis(p8, pi / 4);
    
    nNode = 9;
    Nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    
    POS = [p1, p2, p3, p4, p5, p6, p7, p8, p9];
    
    IEN_s = [1, 1, 1, 1, 1, 1, 1, 1;
             2, 3, 4, 5, 6, 7, 8, 9;
             3, 4, 5, 6, 7, 8, 9, 2];
    
    ID = [1, 0, 0, 0, 0, 0, 0, 0, 0;
          2, 0, 0, 0, 0, 0, 0, 0, 0;
          3, 0, 0, 0, 0, 0, 0, 0, 0];
    
    LM = get_LM(ID, IEN_s, Nodes);
    
    Diri_normal = zeros(3, 9);
    Diri_mcn = zeros(3, 9);
    for nn = 2 : 9
        point = POS(:, nn);
        normal = point / norm(point);
        Diri_normal(:, nn) = normal;
        mcn = normal / 2;
        Diri_mcn(:, nn) = mcn;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_eq = 3;

M = zeros(n_eq, n_eq);
Sx = zeros(n_eq, 1);
b = zeros(n_eq, 1);

hh = 0;
for ee = 1 : n_Elem
    node_ele = zeros(3, n_EN);
    for jj = 1 : n_EN
        node_ele(:, jj) = POS(:, IEN_s(jj, ee));
    end

    triangle = node_ele(:, 1:3);
    hh_ele = circumcircle_3D(triangle);
    if hh_ele > hh
        hh = hh_ele;
    end

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

        G_FFF = F_hat' * F_hat;

        for aa = 1 : n_EN           
            Na = TriBasis(Elem_degree, aa, 0, 0, tqp(:, qua));
            Phi_matrix(1 : 3, 3 * aa - 2 : 3 * aa) = [Na,  0,  0;
                                                       0, Na,  0;
                                                       0,  0, Na];

            % (For mean curvature vector)
            Na_vec = [TriBasis(Elem_degree, aa, 1, 0, tqp(:,qua));
                     TriBasis(Elem_degree, aa, 0, 1, tqp(:,qua))];
            for mm = 1 : 3
                vec1 = zeros(3, 1);
                vec2 = zeros(3, 1);
                Xm_vec = F_hat(mm, :)';
                for nn = 1 : 3
                    Xn_vec = F_hat(nn, :);
                    vec1(nn, 1) = Xn_vec / G_FFF * Na_vec;
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
            Sx(LM_a) = Sx(LM_a) + gamma1 * Sx_ele(aa);
            for bb = 1 : 3 * n_EN
                LM_b = LM(bb, ee);
                if LM_b > 0
                    M(LM_a, LM_b) = M(LM_a, LM_b) + gamma1 * M_ele(aa, bb);
                else
                    % Dirichlet BC
                    bb_node = floor((bb-1) / 3) + 1;
                    g_normal = Diri_normal(:, IEN_s(bb_node, ee));
                    g_mcn = Diri_mcn(:, IEN_s(bb_node, ee));

                    bb_direc = mod(bb, 3);
                    if bb_direc == 0
                        bb_direc = 3;
                    end
                    b(LM_a) = b(LM_a) - gamma1 * M_ele(aa, bb) * g_normal(bb_direc);
                    Sx(LM_a) = Sx(LM_a) - gamma1 * M_ele(aa, bb) * g_mcn(bb_direc);
                end
            end
        end
    end
end

nh_normal = M \ b

nh_mcn = M \ Sx
