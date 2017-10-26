function [idx dist] = rangesearch(P, Q, r)

r2 = r*r;

rowsP = size(P,1);
idx = cell(rowsP,1);
dist = cell(rowsP,1);
for i=1:rowsP
    p = P(i,:);
    [idxp distp] = range(p, Q, r2);
    idx{i} = idxp;
    dist{i} = sqrt(distp);   
end

function [idx dist] = range(p, Q, rsq)
idx = [];
dist = [];
for i = 1:size(Q,1)
    q = Q(i,:);
    diff = q-p;
    distsq = diff*diff';
    if distsq < rsq
        idx = [idx i];
        dist = [dist distsq];
    end
    [dist idxperm] = sort(dist);
    idx = idx(idxperm);
end