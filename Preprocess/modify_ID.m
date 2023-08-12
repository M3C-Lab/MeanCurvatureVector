function new_ID = modify_ID(old_ID, Diri_node)
% To modify the ID array according to the info of Dirichlet nodes.

nNode = size(old_ID, 2);
nDiri = length(Diri_node);

for nn = 1 : nDiri
    old_ID(:, Diri_node(nn)) = [0; 0; 0];
    
    for ii = (Diri_node(nn) + 1) : nNode
        if old_ID(1, ii) ~= 0
            old_ID(:, ii) = old_ID(:, ii) - [3; 3; 3];           
        end
    end
end
new_ID = old_ID;
end

