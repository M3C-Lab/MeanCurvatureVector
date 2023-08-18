clc;clear all;close all;

% Consider a problem with L-B operator defined on R3 space: -(L-B oprator)f = u,
% but we focus on a small domain that on a plane z = 1.

f = @(x, y, z) (x^5 + y^5 + z^5) / 20;
f_x = @(x, y, z) [x^4; y^4; z^4] / 4;
exact_u = @(x, y, z) -(x^3 + y^3);

Elem_degree = 1;

n_EN = 3;

[tqp, wtqp, ntqp] = TriQuad(6);

n_Elem = 5;

phi = 2 * pi / n_Elem;
Radius = 0.001;

nNode = n_Elem + 1;
POS = zeros(3, nNode);
ID = zeros(1, nNode);
ID(1, 1) = 1;

IEN_s = zeros(3, n_Elem);
IEN_s(:, 1) = [1; 2; 3];

LM = zeros(3, n_Elem);

Diri_g = zeros(1, nNode);

for nn = 2 : nNode
    POS(:, nn) = [Radius * cos((nn - 2) * phi) + 1; Radius * sin((nn - 2) * phi); 1];
    Diri_g(1, nn) = exact_u(POS(1, nn), POS(2, nn), POS(3, nn));
end
POS(1, 1) = POS(1, 1) + 1;
POS(2, 1) = POS(2, 1);
POS(3, 1) = 1;

for ee = 1 : n_Elem
    IEN_s(:, ee) = [1; ee + 1; ee + 2];
    LM(1, ee) = 1;
end
IEN_s(3, n_Elem) = 2;

n_eq = 1;

M = zeros(n_eq, n_eq);
Sf = zeros(n_eq, 1);

for ee = 1 : n_Elem
    node_ele = zeros(3, n_EN);
    f_ele = zeros(n_EN, 1);
    f_x_ele = zeros(3, n_EN);
    for aa = 1 : n_EN
        node_ele(:, aa) = POS(:, IEN_s(aa, ee));
        f_ele(aa, 1) = f(node_ele(1, aa), node_ele(2, aa), node_ele(3, aa));
        f_x_ele(:, aa) = f_x(node_ele(1, aa), node_ele(2, aa), node_ele(3, aa));
    end

    M_ele = zeros(n_EN, n_EN);
    Sf_ele = zeros(n_EN, 1);

    for qua = 1 : ntqp
        Phi_matrix = zeros(1, n_EN); % It is actually a vector;

        Sf_qua = zeros(n_EN, 1);

        F_hat = zeros(3, 2);

        N_xi = zeros(2, 3);

        for ii = 1 : 3
            F_hat(ii, :) = [Grad(Elem_degree, ii, 1, node_ele, tqp(:, qua)),...
                            Grad(Elem_degree, ii, 2, node_ele, tqp(:, qua))];

            N_xi(:, ii) = [TriBasis(Elem_degree, ii, 1, 0, tqp(:, qua));
                TriBasis(Elem_degree, ii, 0, 1, tqp(:, qua))];
        end
        G_FFF = F_hat' * F_hat;
        
        xyz_qua = node_ele * tqp(:, qua);

        f_xi_qua = F_hat' * f_x(xyz_qua(1), xyz_qua(2), xyz_qua(3));

        for aa = 1 : n_EN
            Phi_matrix(1, aa) = TriBasis(Elem_degree, aa, 0, 0, tqp(:, qua));

            Na_xi = [TriBasis(Elem_degree, aa, 1, 0, tqp(:, qua));
                TriBasis(Elem_degree, aa, 0, 1, tqp(:, qua))];

            vec1 = zeros(3, 1);
            vec2 = zeros(3, 1);
            for nn = 1 : 3
                Xn_xi = F_hat(nn, :);
                vec1(nn, 1) = Xn_xi / G_FFF * Na_xi;
                vec2(nn, 1) = Xn_xi / G_FFF * f_xi_qua;
            end

            Sf_qua(aa, 1) = dot(vec1, vec2);
        end

        Phi_x_matrix = F_hat / G_FFF * N_xi;

        [normal_f, Jacobian] = get_normal(Elem_degree, node_ele, tqp(:, qua));

        M_ele = M_ele + Jacobian * wtqp(qua) * (Phi_matrix' * Phi_matrix);

        % Sf_ele = Sf_ele + Jacobian * wtqp(qua) * Sf_qua;
        % Approach 1: Real 1st order derivative.
        
        Sf_ele = Sf_ele + Jacobian * wtqp(qua) * (Phi_x_matrix' * Phi_x_matrix) * f_ele;
        % Approach 2: Approximated 1st order derivative with value of f.

        %f_x_qua = f_x_ele * tqp(:, qua);
        %Sf_ele = Sf_ele + Jacobian * wtqp(qua) * Phi_x_matrix' * f_x_qua;
        % Approach 3: Another approximated derivative with value of f_x.
    end

    for aa = 1 : n_EN
        LM_a = LM(aa, ee);
        if LM_a > 0
            Sf(LM_a) = Sf(LM_a) + Sf_ele(aa);
            for bb = 1 : n_EN
                LM_b = LM(bb, ee);
                if LM_b > 0
                    M(LM_a, LM_b) = M(LM_a, LM_b) + M_ele(aa, bb);
                else
                    g = Diri_g(1, IEN_s(bb, ee));
                    Sf(LM_a) = Sf(LM_a) - M_ele(aa, bb) * g;
                end
            end
        end
    end
end

uh = M \ Sf