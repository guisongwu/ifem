function Hdm = stokes_hessian(node, elem, stokes_info, bdFlag, m, um, us1, m1, option);

  
    global slope;

    nx = stokes_info.nx;
    ny = stokes_info.ny;
    n1 = nx - 1;
    xtop = stokes_info.xtop;
    xbot = stokes_info.xbot;
    IUxBot = stokes_info.IUxBot;
    IUyBot = stokes_info.IUyBot;
    IUxTop = stokes_info.IUxTop;

    % size(m)
    % size(um)
    % size(us1)
    % size(m1)
    if size(m1, 1) == n1
        m1 = extend_mid(m1);
    end
    
    assert(n1*2 == length(m));
    Hdm = zeros(n1, 1);
    EI = eye(n1);

    % if size(m, 1) == n1
    %     m = extend_mid(m);
    % end

    
    % 
    % Incremental: 
    %    (dLdu) du = - dLdm m1 
    % 
    pde_du.f = 0;
    pde_du.fp = 0;
    pde_du.g_N = [xtop*0, xtop*0];
    pde_du.g_R = linearize_bot(m);
    %pde_du.g_RN = [-linearize_bot(dLdm1(:,1)), -linearize_bot(dLdm1(:,2))];
    pde_du.g_RN = [-um(IUxBot), m1, -um(IUyBot), m1]; 
    pde_du.g_Dn = [xbot*0, xbot*0];
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_du,option);
    % [du1] = poisson_robin2D_mpfr(x, Z(:)*0, 0*x, -dLdm1, m);
    du1 = soln.u;

    
    % 
    % 
    % Gauss Newton 
    % 
    % 
    % [us3] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
    %                              0*x, m);
    %pde_adj2 = pde;
    pde_adj3.f = 0;
    pde_adj3.fp = 0;
    pde_adj3.g_N = - [2 * linearize_top(du1(IUxTop)), xtop*0];
    pde_adj3.g_R = linearize_bot(m);
    pde_adj3.g_RN = [xbot*0, xbot*0];
    pde_adj3.g_Dn = [xbot*0, xbot*0];
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj3,option);
    us3 = soln.u;

    for i = 1:n1
        m2 = extend_mid(EI(:,i));

        % main term
        term = integral_robin_P2(node, elem, bdFlag, ...
                                 [us3(IUxBot), us3(IUyBot)], ...
                                 [um(IUxBot), um(IUyBot)], ...
                                 [m2, m2], option);

        Hdm(i) = term;
    end

    Hdm = Hdm + stokes_info.gamma_stab * stokes_info.Mstab * m1(1:n1);
    % size(Hdm)
    
end

