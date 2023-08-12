function IEN_v = get_IEN_v(msh, Elem_degree)
% To get the IEN array of volumetric elements.

if Elem_degree == 1
    IEN_v = zeros(4, msh.nbTets);
    temp = 1;
    for jj = msh.nbLines + msh.nbTriangles + 1 : msh.nbElm
        for kk = 1 : 4
            IEN_v(kk, temp) = msh.ELE_NODES(jj, kk);
        end
        temp = temp + 1;
    end
    
elseif Elem_degree == 2
    IEN_v = zeros(10, msh.nbTets10);
    temp = 1;
    for jj = msh.nbLines3 + msh.nbTriangles6 + 1 : msh.nbElm
        for kk = 1 : 10
            IEN_v(kk, temp) = msh.ELE_NODES(jj, kk);
        end
        temp = temp + 1;
    end
end

end

% EOF