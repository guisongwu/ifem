function [int_val] = integral_robin_P2(node,elem,bdFlag,func1,func2,func3,option)

    if ischar(func2)
        assert(strcmp(func2, 'repeat'));
        func2 = func1;
    end
    
    N = size(node,1);
    [elem2dof,edge,bdDof] = dofP2(elem);
    elem2edge = elem2dof(:,4:6)-N;
    
    isRobin(elem2edge(bdFlag(:)==3)) = true;
    Robin     = edge(isRobin,:);
    RobinIdx = find(isRobin);        
    Robin   = edge(isRobin,:);
    nRobin = size(Robin, 1);

    % if ~isfield(option,'quadorder')
    %     option.quadorder = 4;   
    % end
    quadorder1d = 6;
    [lambda,w] = quadpts1(quadorder1d);
    nQuad = size(lambda,1);

    % quadratic bases (1---3---2)
    bdphi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
    bdphi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
    bdphi(:,3) = 4*lambda(:,1).*lambda(:,2);

    % length of edge
    ve = node(Robin(:,1),:) - node(Robin(:,2),:);
    edgeLength = sqrt(sum(ve.^2,2));
    
    assert(~isempty(func1));
    assert(~isempty(func2));
    % assert(iscolumn(func1));
    % assert(iscolumn(func2));
    assert(size(func1, 2) == 1 || size(func1, 2) == 2);
    assert(size(func2, 2) == 1 || size(func2, 2) == 2);
    ncomp = size(func1, 2);

    ge = zeros(size(Robin,1), 1);
    
    for pp = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(pp,1)*node(Robin(:,1),:)+lambda(pp,2)*node(Robin(:,2),:);

        f1 = func1([1:nRobin], :);
        f2 = func1([2:nRobin, 1], :);
        f3 = func1([nRobin+1:2*nRobin], :);

        g1 = func2([1:nRobin], :);
        g2 = func2([2:nRobin, 1], :);
        g3 = func2([nRobin+1:2*nRobin], :);

        if ~isempty(func3)
            h1 = func3([1:nRobin], :);
            h2 = func3([2:nRobin, 1], :);
            h3 = func3([nRobin+1:2*nRobin], :);
        end
        
        val = zeros(size(Robin,1), 1);
        
        for k = 1:ncomp
            vf = f1(:,k) .*  bdphi(pp,1) ...
                 + f2(:,k) .*  bdphi(pp,2) ...
                 + f3(:,k) .*  bdphi(pp,3);

            vg = g1(:,k) .*  bdphi(pp,1) ...
                 + g2(:,k) .*  bdphi(pp,2) ...
                 + g3(:,k) .*  bdphi(pp,3);

            if isempty(func3)
                val = val + vf .* vg;
            else
                vh = h1(:,k) .*  bdphi(pp,1) ...
                     + h2(:,k) .*  bdphi(pp,2) ...
                     + h3(:,k) .*  bdphi(pp,3);
                val = val + vf .* vg .* vh;
            end
        end

        ge(:) = ge(:) + w(pp)*edgeLength.*val(:);
    end

    int_val= sum(ge);
end
