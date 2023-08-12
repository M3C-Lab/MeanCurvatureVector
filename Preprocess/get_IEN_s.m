function IEN_s = get_IEN_s(msh, Elem_degree, surf_name)
% To get the IEN array of surface elements.

for I = 1 : size(msh.PhyGrp, 1)
    if strcmp(surf_name, msh.PhyGrp{I, 3})
        PG_number = msh.PhyGrp{I, 2};
        break;
    end
end

nElem_s = 0;

for ii = 1 : msh.nbElm
    if msh.ELE_TAGS(ii, 1) == PG_number
        nElem_s = nElem_s + 1;
    end
end

temp = 1;
if Elem_degree == 1
    IEN_s = zeros(3, nElem_s);
    for ii = 1 : msh.nbTriangles
        if msh.TRIANGLES(ii, 4) == PG_number
            for jj = 1 : 3
                IEN_s(jj, temp) = msh.TRIANGLES(ii, jj);
            end
            temp = temp + 1;
        end
    end
    
elseif Elem_degree == 2
    IEN_s = zeros(6, nElem_s);
    for ii = 1 : msh.nbTriangles6
        if msh.TRIANGLES6(ii, 7) == PG_number
            for jj = 1 : 6
                IEN_s(jj, temp) = msh.TRIANGLES6(ii, jj);
            end
            temp = temp + 1;
        end
    end
    
end


end

% EOF

