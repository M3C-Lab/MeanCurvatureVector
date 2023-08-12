function [number, new_qp] = find_sidenode_number(side, IEN_ee, qp)
% Input:
%   side: side_node column vector.
%   IEN_ee: IEN_s column vecter.
%   qp: quadrature point on line element between (0, 1).
% Output:
%   number: the node number in the element of that side.
%   new_qp: quadrature point of area coordinate on the same point in the element.

number1 = find(side(1) == IEN_ee);
number2 = find(side(2) == IEN_ee);

number = [number1; number2];

new_qp = zeros(3, length(qp));
for qua = 1 : length(qp)
    new_qp(number1, qua) = qp(qua);
    new_qp(number2, qua) = 1 - qp(qua);
end

end

% EOF

