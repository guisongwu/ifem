function [int_val] = integral_neumann(node,elem,bdFlag,u1,u2,option)

    Neumann = [];
    isNeumann = (bdFlag(:) == 2);
    if any(isNeumann)
        allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        Neumann = allEdge(isNeumann,:);
    end

    Neumann = sort(Neumann, 2);
    Neumann = sortrows(Neumann);

    el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
    if ~isfield(option,'gNquadorder')
        option.gNquadorder = 2;   % default order exact for linear gN
    end

    assert(~isempty(u1));
    assert(~isempty(u2));
    assert(iscolumn(u1));
    assert(iscolumn(u2));

    u1_shift = [u1(2:end); u1(1)];
    u1_elem = [u1 u1_shift];

    u2_shift = [u2(2:end); u2(1)];
    u2_elem = [u2 u2_shift];

    [lambdagN,weightgN] = quadpts1(option.gNquadorder);
    phigN = lambdagN;                 % linear bases
    nQuadgN = size(lambdagN,1);
    ge = zeros(size(Neumann,1), 1);


    for pp = 1:nQuadgN
        % quadrature points in the x-y coordinate
        ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
               + lambdagN(pp,2)*node(Neumann(:,2),:);

        u1_val = lambdagN(pp,1) * u1_elem(:,1) ...
                 + lambdagN(pp,2) * u1_elem(:,2);

        u2_val = lambdagN(pp,1) * u2_elem(:,1) ...
                 + lambdagN(pp,2) * u2_elem(:,2);

        ge(:) = ge(:) + weightgN(pp)*u1_val.*u2_val;
    end

    ge = ge.*el;
    int_val= sum(ge);
    
end
