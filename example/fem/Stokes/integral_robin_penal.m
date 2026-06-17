function [int_val] = integral_robin_penal(node,elem,bdFlag,func1,func2,option)

    
    N = size(node,1);
    [elem2dof,edge,bdDof] = dofP2(elem);
    elem2edge = elem2dof(:,4:6)-N;

    isRobin = false(size(edge, 1), 1);
    isRobin(elem2edge(bdFlag(:)==3)) = true;
    Robin     = edge(isRobin,:);
    RobinIdx = find(isRobin);        
    Robin   = edge(isRobin,:);
    nRobin = size(Robin, 1);

    RobinElem = zeros(size(Robin, 1), 3);
    for i = 1:size(Robin, 1)
        iE = find(elem2edge(:,1) == RobinIdx(i) | elem2edge(:,2) == RobinIdx(i) | elem2edge(:,1) == RobinIdx(i));
        RobinElem(i, 1) = iE;
        RobinElem(i, 2) = find(Robin(i, 1) == elem(iE, :));
        RobinElem(i, 3) = find(Robin(i, 2) == elem(iE, :));
    end


    % if ~isfield(option,'quadorder')
    %     option.quadorder = 4;   
    % end
    quadorder1d = 6;
    [lambda,w] = quadpts1(quadorder1d);
    nQuad = size(lambda,1);

    [DlambdaAll,area] = gradbasis(node,elem); % Dlambda: [Nelem][x|y][lambda[123]]
    Dlambda = zeros(size(Robin, 1), 2, 2);
    for i = 1:size(Robin, 1)
        Dlambda(i, :, 1) = DlambdaAll(RobinElem(i,1), :, RobinElem(i, 2));
        Dlambda(i, :, 2) = DlambdaAll(RobinElem(i,1), :, RobinElem(i, 3));
    end
    
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
    assert(size(func1, 2) == 1 );
    assert(size(func2, 2) == 1 );


    
    ge = zeros(size(Robin,1), 1);
    
    for pp = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(pp,1)*node(Robin(:,1),:)+lambda(pp,2)*node(Robin(:,2),:);

        Dphip(:,:,3) = 4*(lambda(pp,1)*Dlambda(:,:,2)+lambda(pp,2)*Dlambda(:,:,1));
        % Dphip(:,:,5) = 4*(lambda(pp,3)*Dlambda(:,:,1)+lambda(pp,1)*Dlambda(:,:,3));
        % Dphip(:,:,4) = 4*(lambda(pp,2)*Dlambda(:,:,3)+lambda(pp,3)*Dlambda(:,:,2));
        Dphip(:,:,1) = (4*lambda(pp,1)-1).*Dlambda(:,:,1);            
        Dphip(:,:,2) = (4*lambda(pp,2)-1).*Dlambda(:,:,2);            
        %Dphip(:,:,3) = (4*lambda(pp,3)-1).*Dlambda(:,:,3);            
        
        f1 = func1([1:nRobin], :);
        f2 = func1([2:nRobin, 1], :);
        f3 = func1([nRobin+1:2*nRobin], :);

        g1 = func2([1:nRobin], :);
        g2 = func2([2:nRobin, 1], :);
        g3 = func2([nRobin+1:2*nRobin], :);
       
        val = zeros(size(Robin,1), 1);
        
        % vf = f1(:,k) .*  bdphi(pp,1) ...
        %      + f2(:,k) .*  bdphi(pp,2) ...
        %      + f3(:,k) .*  bdphi(pp,3);

        % vg = g1(:,k) .*  bdphi(pp,1) ...
        %      + g2(:,k) .*  bdphi(pp,2) ...
        %      + g3(:,k) .*  bdphi(pp,3);

        vf = f1(:) .*  Dphip(:,:,1) ...
             + f2(:) .*  Dphip(:,:,2) ...
             + f3(:) .*  Dphip(:,:,3);
        % [nRobin][d[x|y]]

        vg = g1(:) .*  Dphip(:,:,1) ...
             + g2(:) .*  Dphip(:,:,2) ...
             + g3(:) .*  Dphip(:,:,3);

        size(vf);
        size(vg);
        
        if 0
            % 
        end
        

        %ge(:) = ge(:) + w(pp)*edgeLength.* (vf(:,1) .* vg(:,1) + vf(:,2) .* vg(:,2))
        ge(:) = ge(:) + w(pp)*edgeLength.* (vf(:,1) .* vg(:,1));
    end

    int_val= sum(ge);
end
