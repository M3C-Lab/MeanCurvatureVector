function facet_normal = get_facet_normal(msh, IEN_v, IEN_s)
% To get facet normal of 1st order triangle elements:
% Output:
%   facet_normal(:, ii): The facet normal of ii-th element.

%%%%% Only for 1st order element! %%%%%

facet_normal = zeros(3, size(IEN_s, 2));

for ii = 1 : size(IEN_s, 2)
    tri = IEN_s(:, ii);
    for jj = 1 : size(IEN_v, 2)
        tet = IEN_v(:, jj);
        tet_nb = 0;
        for kk = 1 : 4
            if tri(1) == tet(kk)
                for mm = 1 : 4
                    if tri(2) == tet(mm)
                        for nn = 1 : 4
                            if tri(3) == tet(nn)                                
                                tet_nb = jj;
                                break;
                            end
                        end
                    end
                    if tet_nb ~= 0
                        break;
                    end
                end
            end
            if tet_nb ~= 0
                break;
            end
        end
        
        
        
        if tet_nb ~= 0
            for kk = 1 : length(tet)
                if tet(kk) == tri(1)
                    tet(kk) = [ ];
                    break;
                end
            end
            for kk = 1 : length(tet)
                if tet(kk) == tri(2)
                    tet(kk) = [ ];
                    break;
                end
            end
            for kk = 1 : length(tet)
                if tet(kk) == tri(3)
                    tet(kk) = [ ];
                    break;
                end
            end
            opposite = tet(1);
            p1 = [msh.POS(IEN_s(1, ii), 1), msh.POS(IEN_s(1, ii), 2), msh.POS(IEN_s(1, ii), 3)]';
            p2 = [msh.POS(IEN_s(2, ii), 1), msh.POS(IEN_s(2, ii), 2), msh.POS(IEN_s(2, ii), 3)]';
            p3 = [msh.POS(IEN_s(3, ii), 1), msh.POS(IEN_s(3, ii), 2), msh.POS(IEN_s(3, ii), 3)]';
            p4 = [msh.POS(opposite, 1), msh.POS(opposite, 2), msh.POS(opposite, 3)]';
            
            edge1 = p2 - p1;
            edge2 = p3 - p1;
            normal = cross(edge1, edge2);
            test_edge = p4 - p1;
            if dot(normal, test_edge) > 0
                normal = -normal;
            end
            normal = normal / norm(normal);
            facet_normal(:, ii) = normal;
            break;
        end
    end
end

end

% EOF

