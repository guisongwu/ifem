    pde = Stokesdata2;
    [node,elem] = squaremesh([0,1,0,1],0.25);
    err = zeros(4,1); h = zeros(4,1);
    for k = 1:4
        [node,elem] = uniformrefine(node,elem);
        [elem2edge,edge] = dofedge(elem);
        uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
        h(k) = 1/sqrt(size(elem,1));
        err(k) = getL2error(node,elem,pde.exactu,uI);
    end
    showrateh(h,err);