function LM = get_LM(ID, IEN, Nodes)
% To get the LM array.

nbENod = size(IEN, 1);
nbEle = size(IEN, 2);

LM = zeros(3 * nbENod, nbEle);
for jj = 1 : nbEle
    for ii = 1 : nbENod
        for kk = 1 : 3
            LM(3*(ii-1)+kk , jj) = ID(kk, find(Nodes == IEN(ii, jj)));
        end
    end
end

end

% EOF