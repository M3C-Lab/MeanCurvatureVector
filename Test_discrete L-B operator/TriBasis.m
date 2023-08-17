function poly = TriBasis(degree, i, der_r, der_s, area_coor)
% To get the basis function or the derivative over a triangular element.
% phys2rst provides the AREA COORDINATES, which are crucial for the basis functions.
% This formulation is on the pages of 167 ~ 168 of Hughes' book.
% Input:
%   degree: the degree of triangular element, 1 ~ 2
%   i: the number of the basis function, 1 ~ 3.
%   der: the derivative required.
%       0 --- the value of basis function.
%       1 --- the 1st partial derivative of the basis function.
%       2 --- the 2nd partial derivative of the basis function.
% Output:
%   poly: The value of basis function or derivative.

r = area_coor(1);
s = area_coor(2);
t = area_coor(3);

if degree == 1
    if i == 1
        if der_r == 0 && der_s == 0
            poly = r;
        elseif der_r == 1 && der_s == 0
            poly = 1;
        elseif der_r == 0 && der_s == 1
            poly = 0;
        else
            poly = 0;
        end

    elseif i == 2
        if der_r == 0 && der_s == 0
            poly = s;
        elseif der_r == 1 && der_s == 0
            poly = 0;
        elseif der_r == 0 && der_s == 1
            poly = 1;
        else
            poly = 0;
        end

    elseif i == 3
        if der_r == 0 && der_s == 0
            poly = t;
        elseif der_r == 1 && der_s == 0
            poly = -1;
        elseif der_r == 0 && der_s == 1
            poly = -1;
        else
            poly = 0;
        end

    else
        disp("TriBasis: Please input appropriate node number.")
    end
    
elseif degree == 2
    if i == 1
        if der_r == 0 && der_s == 0
            poly = r * (2 * r - 1);
        elseif der_r == 1 && der_s == 0
            poly = (4 * r - 1);
        elseif der_r == 0 && der_s == 1
            poly = 0;
        elseif der_r == 2 && der_s == 0
            poly = 4;
        elseif der_r == 0 && der_s == 2
            poly = 0;
        elseif der_r == 1 && der_s == 1
            poly = 0;
        else
            poly = 0;
        end

    elseif i == 2
        if der_r == 0 && der_s == 0
            poly = s * (2 * s - 1);
        elseif der_r == 1 && der_s == 0
            poly = 0;
        elseif der_r == 0 && der_s == 1
            poly = (4 * s - 1);
        elseif der_r == 2 && der_s == 0
            poly = 0;
        elseif der_r == 0 && der_s == 2
            poly = 4;
        elseif der_r == 1 && der_s == 1
            poly = 0;
        else
            poly = 0;
        end

    elseif i == 3
        if der_r == 0 && der_s == 0
            poly = t * (2 * t - 1);
        elseif der_r == 1 && der_s == 0
            poly = 4 * r + 4 * s - 3;
        elseif der_r == 0 && der_s == 1
            poly = 4 * r + 4 * s - 3;
        elseif der_r == 2 && der_s == 0
            poly = 4;
        elseif der_r == 0 && der_s == 2
            poly = 4;
        elseif der_r == 1 && der_s == 1
            poly = 4;
        else
            poly = 0;
        end
        
    elseif i == 4
        if der_r == 0 && der_s == 0
            poly = 4 * r * s;
        elseif der_r == 1 && der_s == 0
            poly = 4 * s;
        elseif der_r == 0 && der_s == 1
            poly = 4 * r;
        elseif der_r == 2 && der_s == 0
            poly = 0;
        elseif der_r == 0 && der_s == 2
            poly = 0;
        elseif der_r == 1 && der_s == 1
            poly = 4;
        else
            poly = 0;
        end
        
    elseif i == 5
        if der_r == 0 && der_s == 0
            poly = 4 * s * t;
        elseif der_r == 1 && der_s == 0
            poly = -4 * s;
        elseif der_r == 0 && der_s == 1
            poly = -8 * s - 4 * r + 4;
        elseif der_r == 2 && der_s == 0
            poly = 0;
        elseif der_r == 0 && der_s == 2
            poly = -8;
        elseif der_r == 1 && der_s == 1
            poly = -4;
        else
            poly = 0;
        end
        
    elseif i == 6
        if der_r == 0 && der_s == 0
            poly = 4 * r * t;
        elseif der_r == 1 && der_s == 0
            poly = -8 * r - 4 * s + 4;
        elseif der_r == 0 && der_s == 1
            poly = -4 * r;
       elseif der_r == 2 && der_s == 0
            poly = -8;
        elseif der_r == 0 && der_s == 2
            poly = 0;
        elseif der_r == 1 && der_s == 1
            poly = -4;
        else
            poly = 0;
        end

    else
        disp("TriBasis: Please input appropriate node number.")
    end
end

end

% EOF
