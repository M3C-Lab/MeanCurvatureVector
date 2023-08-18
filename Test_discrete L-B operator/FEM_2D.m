clc; clear all; close all;

f = @(x, y) x^3 + y^3;
f_x = @(x, y) [3*x^2 ; 3 * y^2];
exact_u = @(x, y) - 6 * x  - 6 * y;

Elem_degree = 1;

n_EN = 3;

[tqp, wtqp, ntqp] = TriQuad(6);

p1 = [0; 0];

n_Elem = 6;

phi = 2 * pi / n_Elem;
Radius = 1;

nNode = n_Elem + 1;
POS = zeros(2, nNode);
ID = zeros(1, nNode);
ID(1, 1) = 1;

IEN_s = zeros(3, n_Elem);
IEN_s(:, 1) = [1; 2; 3];

LM = zeros(3, n_Elem);

Diri_g = zeros(1, nNode);

for nn = 2 : nNode
    POS(:, nn) = [Radius * cos((nn - 2) * phi) + 1; Radius * sin((nn - 2) * phi) + 0.25];
    Diri_g(1, nn) = exact_u(POS(1, nn), POS(2, nn));
end
POS(1, 1) = POS(1, 1) + 1;
POS(2, 1) = POS(2, 1) + 0.25;

for ee = 1 : n_Elem
    IEN_s(:, ee) = [1; ee + 1; ee + 2];
    LM(1, ee) = 1;
end
IEN_s(3, n_Elem) = 2;

n_eq = 1;

M = zeros(n_eq, n_eq);
Sf = zeros(n_eq, 1);

for ee = 1 : n_Elem
    node_ele = zeros(2, n_EN);
    f_ele = zeros(3, n_EN);
    f_x_ele = zeros(2, n_EN);
    for aa = 1 : n_EN
        node_ele(:, aa) = POS(:, IEN_s(aa, ee));
        f_ele(aa, 1) = f(node_ele(1, aa), node_ele(2, aa));
        f_x_ele(:, aa) = f_x(node_ele(1, aa), node_ele(2, aa));
    end

    M_ele = zeros(n_EN, n_EN);
    Sf_ele = zeros(n_EN, 1);

    for qua = 1 : ntqp
        F_hat = zeros(2, 2);
        
        for ii = 1 : 2
            F_hat(ii, :) = [Grad(Elem_degree, ii, 1, node_ele, tqp(:, qua)),...
                            Grad(Elem_degree, ii, 2, node_ele, tqp(:, qua))];
        end

        Jacobian = det(F_hat);

        xy_qua = node_ele * tqp(:, qua);

        f_x_qua = f_x(xy_qua(1), xy_qua(2));

        Phi_matrix = zeros(1, n_EN);

        N_xi = zeros(2, 3);

        for aa = 1 : n_EN
            Phi_matrix(1, aa) = TriBasis(Elem_degree, aa, 0, 0, tqp(:, qua));

            N_xi(:, aa) = [TriBasis(Elem_degree, aa, 1, 0, tqp(:, qua));
                TriBasis(Elem_degree, aa, 0, 1, tqp(:, qua))];
        end

        Phi_x_matrix = inv(F_hat') *  N_xi;

        M_ele = M_ele + Jacobian * wtqp(qua) * (Phi_matrix' * Phi_matrix);

        % Sf_ele = Sf_ele + Jacobian * wtqp(qua) * (Phi_x_matrix' * f_x_qua);
        % Real 1st order derivative.
    
        %Sf_ele = Sf_ele + Jacobian * wtqp(qua) * (Phi_x_matrix' * Phi_x_matrix) * f_ele;
        % Approximated 1st order derivative.

        f_x_qua = f_x_ele * tqp(:, qua);
        Sf_ele = Sf_ele + Jacobian * wtqp(qua) * Phi_x_matrix' * f_x_qua;
        % Approximated 2.
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