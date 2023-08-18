clc; clear all;close all;

% To solve - d2f / dx2 = u with 1st order element.

f = @(x) sin(x);
f_x = @(x) cos(x);
exact_u = @(x) sin(x);

omega_l = 0.0;
omega_r = 10.0;

pp = 1;

n_Elem = 30;

nqp = 6;
[qp, wqp] = Gauss(nqp, -1, 1);

n_np = n_Elem * pp + 1;
n_en = pp + 1;
n_eq = n_np - 2;

IEN = zeros(n_en, n_Elem);
for ee = 1 : n_Elem
    for aa = 1 : n_en
        IEN(aa, ee) = (ee - 1) * pp + aa;
    end
end

hh = (omega_r - omega_l) / n_Elem;

x_coor = omega_l : (hh / pp) : omega_r;

ID = 1 : n_np;
ID(1) = 0;
ID(end) = 0;
for nn = 2 : n_np - 1
    ID(nn) = ID(nn) - 1;
end

uh = [exact_u(omega_l); zeros(n_eq, 1); exact_u(omega_r)];

M = zeros(n_eq, n_eq);
Sf = zeros(n_eq, 1);

for ee = 1 : n_Elem
    M_ele = zeros(n_en, n_en);
    Sf_ele = zeros(n_en, 1);

    x_ele = zeros(n_en, 1);
    f_ele = zeros(n_en, 1);
    f_x_ele = zeros(n_en, 1);

    for aa = 1 : n_en
        x_ele(aa) = x_coor(IEN(aa, ee));
        f_ele(aa) = f(x_ele(aa));
        f_x_ele(aa) = f_x(x_ele(aa));
    end

    for qua = 1 : nqp
        x_qua = 0.0;
        dx_dxi = 0.0;
        
        for aa = 1 : n_en
            x_qua = x_qua + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end

        dxi_dx = 1.0 / dx_dxi;

        Phi_matrix = zeros(1, n_en);
        Phi_x_matrix = zeros(1, n_en);

        for aa = 1 : n_en
            Phi_matrix(1, aa) = PolyBasis(pp, aa, 0, qp(qua));
            Phi_x_matrix(1, aa) = PolyBasis(pp, aa, 1, qp(qua)) * dxi_dx;
        end

        M_ele = M_ele + wqp(qua) * dx_dxi * (Phi_matrix' * Phi_matrix);
        Sf_ele = Sf_ele + wqp(qua) * dx_dxi * (Phi_x_matrix' * Phi_x_matrix) * f_ele;
    end

    for aa = 1 : n_en
        LM_a = ID(IEN(aa, ee));
        if LM_a > 0
            Sf(LM_a) = Sf(LM_a) + Sf_ele(aa);
            for bb = 1 : n_en
                LM_b = ID(IEN(bb, ee));
                if LM_b > 0
                    M(LM_a, LM_b) = M(LM_a, LM_b) + M_ele(aa, bb);
                else
                    x_b = x_coor(IEN(bb, ee));
                    Sf(LM_a) = Sf(LM_a) - M_ele(aa, bb) * exact_u(x_b);
                end
            end
        end
    end
end

uh_solution = M \ Sf;
uh(2 : n_np - 1) = uh_solution;

figure(1);
plot(x_coor, exact_u(x_coor), '-b', 'LineWidth', 2);
hold on;
plot(x_coor, uh, '*r', 'LineWidth', 2);
xlabel('x', 'FontSize', 16);
ylabel('u(x)','FontSize', 16);
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 16;
lgd = legend('exact solution','u_h');
lgd.FontSize = 16;
