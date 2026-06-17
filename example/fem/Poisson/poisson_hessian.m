function Hdm = poisson_hessian(node, elem, bdFlag, m, Um, Us1, m1, option);

    global slope;

    [ny, nx] = size(Um);
    % xx = unique(node(:, 1));
    % yy = unique(node(:, 2) + node(:,1) * slope);
    % nx = length(xx);
    % ny = length(yy);
    n1 = nx - 1;
    %h = xx(2) - xx(1);

    X = reshape(node(:,1), ny, nx);
    Y = reshape(node(:,2), ny, nx);
    xx = X(1, 1:nx)';

    xbot = xx(1:n1);
    xtop = xbot;

    assert(n1 == length(m));
    Hdm = zeros(n1, 1);

    % assert(size(Um, 1) == ny);
    % assert(size(Um, 2) == nx);
    
    % 
    % Incremental: 
    %    (dLdu) du = - dLdm m1 
    % 
    dLdm1 = Um(1, 1:n1)' .* m1;
    %pde_du = pde;
    pde_du.f = 0;
    pde_du.g_N = xtop*0;
    pde_du.g_R = m;
    pde_du.g_RN = -dLdm1;
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_du,option);
    % [du1] = poisson_robin2D_mpfr(x, Z(:)*0, 0*x, -dLdm1, m);
    du1 = soln.u;
    Du1 = reshape(du1, ny, nx);
    EI = eye(n1);

    
    if option.use_newton
        
        % 
        % Incremental adj: 
        %    L^t u^s = - dXidu^t 
        % 
        % [us2] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
        %                              Us1(1,:) .* m1', m);
        %pde_adj2 = pde;
        pde_adj2.f = 0;
        pde_adj2.g_N = - 2 * Du1(ny, 1:n1)';
        pde_adj2.g_R = m;
        pde_adj2.g_RN = - Us1(1,1:n1)' .* m1; 
        [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj2,option);
        us2 = soln.u;
        Us2 = reshape(us2, ny, nx);
        
        for i = 1:n1
            m2 = EI(:,i);
            
            % main term
            dLdm2 = Um(1, 1:n1)' .* m2;        
            % term1 = dot(Us2(1,1:n1)', dLdm2);
            term1 = integral_robin(node, elem, bdFlag, Us2(1,1:n1)', dLdm2,option);
            

            % second term is much smaller
            dLdm2_du1 = Du1(1, 1:n1)' .* m2;
            % term2 = dot(Us1(1,1:n1)', dLdm2_du1);
            term2 = integral_robin(node, elem, bdFlag, Us1(1,1:n1)', dLdm2_du1, option);

            term3 = 0;
            term = term1 + term2 + term3;

            % if i == 1 && m1(1) == 1
            %     fprintf('%e\n', [norm(Us1)
            %                      norm(Um)
            %                      norm(Du1)
            %                      norm(dLdm2)
            %                      norm(Us2)]);
            %     [term1 term2 term]
            % end

            Hdm(i) = term;
        end


    else

        % 
        % 
        % Gauss Newton 
        % 
        % 
        % [us3] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
        %                              0*x, m);
        %pde_adj2 = pde;
        pde_adj2.f = 0;
        pde_adj2.g_N = - 2 * Du1(ny, 1:n1)';
        pde_adj2.g_R = m;
        pde_adj2.g_RN = xbot*0;
        [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj2,option);
        us3 = soln.u;
        Us3 = reshape(us3, ny, nx);

        for i = 1:n1
            m2 = EI(:,i);

            % main term
            dLdm2 = Um(1, 1:n1)' .* m2;        
            term = integral_robin(node, elem, bdFlag, Us3(1,1:n1)', dLdm2, option);

            Hdm(i) = term;
        end

    end

    Hdm = Hdm + 1e-9 * m1;
    
end

