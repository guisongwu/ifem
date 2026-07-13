function [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option)
%% STOKESP2P1 Stokes equation: P2-P1 Taylor-Hood elements.
%
%   [soln,eqn,info] = STOKESP2P1(node,elem,pde,bdFlag) use quadratic and piceswise
%   linear elements to approximate velocity u and pressure p, repectively.
% 
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%   with 
%       Dirichlet boundary condition        u = g_D  on \Gamma_D, 
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.
%
%   It is a choice of option.fem in Stokes. Please read Stokes for more
%   information on the input and output.
%
% See also Stokes, Poisson, StokesP2P1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

    global slope dbg_on;

    if ~exist('option','var'), option = []; end

    if option.use_slip
        use_slip = true;
        if option.verb > 0 
            warning('Use slip');
        end
    else
        use_slip = false;
    end

    t = cputime;
    %% Construct Data Structure
    [elem2dof,edge,bdDof] = dofP2(elem);
    N = size(node,1);  NT = size(elem,1);  Nu = N+size(edge,1);   Np = N;


    periodic = option.periodic;
    if periodic
        if option.verb >= 1
            fprintf(2, 'Waring: force periodic\n');
        end
        unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];

        % 
        % u
        % 
        I0 = find(unode(:,1) == 0);
        I0 = sort(I0);
        I1 = find(unode(:,1) == 1);
        I1 = sort(I1);
        assert(norm(unode(I0, 2) - (unode(I1, 2) + slope*(unode(I1, 1))) ) < 1e-12);
        
        Iu_period = [1:2*Nu]';
        Iu_period(I1) = I0;
        Iu_period(Nu+I1) = Nu+I0;
        Mu_period = [I0 I1; I0+Nu I1+Nu];
        Nu_period = length(I1);
        if option.verb >= 1
            Iu_map = [I0 I1];
        end

        %  
        % p
        % 
        I0 = find(node(:,1) == 0);
        I0 = sort(I0);
        I1 = find(node(:,1) == 1);
        I1 = sort(I1);
        assert(norm(node(I0, 2) - (node(I1, 2) + slope*(node(I1, 1))) ) < 1e-12);
        
        Ip_period = [1:Np]';
        Ip_period(I1) = I0;
        Np_period = length(I1);
        Mp_period = [I0 I1];
        if option.verb >= 1
            Ip_map = [I0 I1];
        end

        I_period = [Iu_period; 2*Nu + Ip_period]; % map global dof to global periodic dof
        M_period = [Mu_period; 2*Nu + Mp_period]; % pairs of periodic dof
        N_period = 2*Nu_period + Np_period;

        Ndof = 2*Nu+Np;


        % check
        % unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
        % allnode = [unode; unode; node];
        % allnode(I_period, 1) - allnode(I_period, 1);
        assert(sum(I_period(M_period(:,2)) == M_period(:,1) - 1) == 0);
        assert(sum(I_period(M_period(:,1)) == M_period(:,1) - 1) == 0);
        
        %error('dbg');
        clear I0 I1 
    end




    
    %% Compute geometric quantities and gradient of local basis
    [Dlambda,area] = gradbasis(node,elem);
    % Dlambda(nelem, d[xy], lambda[123]): dlambda / dx

    if option.verb > 0
        fprintf('\nStokes periodic\n');
        fprintf('node %d\n', N);
        fprintf('elem %d\n', NT);
        fprintf('Nu   %d\n', Nu);
        fprintf('Np   %d\n', N);
        fprintf('E2D  %d\n', size(elem2dof));
        %fprintf('Dlambda   %d\n', size(Dlambda));
    end


    %% Assemble stiffness matrix for Laplace operator
    % generate sparse pattern
    ii = zeros(6*6*NT,1); jj = zeros(6*6*NT,1); 
    index = 0;
    for i = 1:6
        for j = 1:6
            ii(index+1:index+NT) = double(elem2dof(:,i)); 
            jj(index+1:index+NT) = double(elem2dof(:,j));  
            index = index + NT;
        end
    end

    % quadrature points
    if ~isfield(pde,'nu'), pde.nu = []; end
    if ~isfield(option,'quadorder')
        % constant viscosity
        option.quadorder = 4;        % default order
        if ~isempty(pde.nu) && isnumeric(pde.nu) % numerical viscosity
            option.quadorder = 4;    % exact for linear diffusion coefficient
        end
    end
    quadorder1d = 6;
    [lambda, w] = quadpts(option.quadorder);
    nQuad = size(lambda,1);
    % compute non-zeros
    sA = zeros(6*6*NT,nQuad);
    for p = 1:nQuad
        % Dphi at quadrature points
        Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
        Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
        Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
        Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
        Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);

        % Dphip[nelem][x|y][nbas]
        % size(Dphip)
        % error('stop')
        
        index = 0;
        for i = 1:6
            for j = 1:6
                Aij = 0;

                if isempty(pde.nu) || isnumeric(pde.nu)
                    Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
                else
                    pxy = lambda(p,1)*node(elem(:,1),:) ...
                          + lambda(p,2)*node(elem(:,2),:) ...
                          + lambda(p,3)*node(elem(:,3),:);
                    Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.d(pxy);
                end
                if ~isempty(pde.nu) && (pde.nu~=1)
                    Aij = pde.nu*Aij;
                end
                
                Aij = Aij.*area;
                sA(index+1:index+NT,p) = Aij;
                index = index + NT;
            end                         % end j
        end                             % end i
    end
    sA = sum(sA,2);                     % quad sum_q

    % assemble the matrix
    
    % diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    % A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Nu,Nu);
    % AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Nu,Nu);
    % A = A + AU + AU';
    A =  sparse(ii,jj,sA,Nu,Nu);
    A = blkdiag(A,A);
    clear Aij ii jj sA


    
    %% Assemble the matrix for divergence operator
    Dx = sparse(Np,Nu);
    Dy = sparse(Np,Nu);
    [lambda, w] = quadpts(2); % (div(P2), P1) is P2
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % Dphi at quadrature points
        Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
        Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
        Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
        Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
        Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));    
        for i = 1:6 
            for j = 1:3
                Dxij = 0;
                Dyij = 0;
                Dxij = Dxij + w(p)*Dphip(:,1,i).*lambda(p,j);
                Dyij = Dyij + w(p)*Dphip(:,2,i).*lambda(p,j);
                Dx = Dx + sparse(elem(:,j),double(elem2dof(:,i)),Dxij.*area,Np,Nu);
                Dy = Dy + sparse(elem(:,j),double(elem2dof(:,i)),Dyij.*area,Np,Nu);
            end
        end
    end
    B = -[Dx Dy];
    clear Dxij Dyij Dx Dy Dphip



    %% Assemble right hand side
    f1 = zeros(Nu,1);
    f2 = zeros(Nu,1);
    if ~isfield(option,'quadorder')
        option.quadorder = 4;   % default order
    end
    if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
        pde.f = [];
        f = zeros(Nu*2,1);
    end
    if ~isempty(pde.f) 
        % quadrature points in the barycentric coordinate
        [lambda,weight] = quadpts(option.quadorder);
        % basis values at quadrature points
        phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
        phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
        phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
        phi(:,4) = 4*lambda(:,2).*lambda(:,3);
        phi(:,5) = 4*lambda(:,3).*lambda(:,1);
        phi(:,6) = 4*lambda(:,1).*lambda(:,2);
        nQuad = size(lambda,1);
        ft1 = zeros(NT,6);
        ft2 = zeros(NT,6);
        for p = 1:nQuad
            % quadrature points in the x-y coordinate
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                  + lambda(p,2)*node(elem(:,2),:) ...
                  + lambda(p,3)*node(elem(:,3),:);
            % function values at quadrature points
            fp = pde.f(pxy);
            % evaluate fp outside.
            for j = 1:6
                ft1(:,j) = ft1(:,j) + fp(:,1).*phi(p,j)*weight(p);
                ft2(:,j) = ft2(:,j) + fp(:,2).*phi(p,j)*weight(p);
            end
        end
        ft1 = ft1.*repmat(area,1,6);
        ft2 = ft2.*repmat(area,1,6);
        f1 = accumarray(elem2dof(:),ft1(:),[Nu 1]);
        f2 = accumarray(elem2dof(:),ft2(:),[Nu 1]);
        f = [f1; f2];
        %norm(f)
    end


    if ~isfield(pde,'fp') || (isreal(pde.fp) && (pde.fp==0))
        pde.fp = [];
        g = zeros(Np,1);
    end
    if ~isempty(pde.fp) 
        % quadrature points in the barycentric coordinate
        [lambda,weight] = quadpts(option.quadorder);

        % basis values at quadrature points
        phi(:,1) = lambda(:,1);
        phi(:,2) = lambda(:,2);
        phi(:,3) = lambda(:,3);
        nQuad = size(lambda,1);
        ft = zeros(NT,3);
        for p = 1:nQuad
            % quadrature points in the x-y coordinate
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                  + lambda(p,2)*node(elem(:,2),:) ...
                  + lambda(p,3)*node(elem(:,3),:);

            % function values at quadrature points
            fpp = -pde.fp(pxy);         % - div(u)

            % evaluate fp outside.
            for j = 1:3
                ft(:,j) = ft(:,j) + fpp(:).*phi(p,j)*weight(p);
            end
        end
        ft = ft.*repmat(area,1,3);
        
        g = accumarray(elem(:),ft(:),[N 1]);
        clear phi ft
    end
    clear phi ft1 ft2 


    %% Boundary Conditions
    [AD,BDt,fD,gD,u,p,ufreeDof,fixedDof,pDof,rot_info] = getbdStokesP2P1;


    %% Record assembeling time
    assembleTime = cputime - t;
    if ~isfield(option,'printlevel'), option.printlevel = 1; end
    if option.printlevel >= 2
        fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
    end


    %% Solve the system of linear equations
    if isempty(ufreeDof), return; end
    if isempty(option) || ~isfield(option,'solver')    % no option.solver
        if length(f)+length(g) <= 1e3  % Direct solver for small size systems
            option.solver = 'direct';
        else          % Multigrid-type  solver for large size systems
            option.solver = 'asmg';
        end
    end
    solver = option.solver;

    % solve
    switch solver
      case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);        
      case 'direct'
        t = cputime;
        % size(AD)
        % size(B)

        if ~use_slip
            bigA = [AD, BDt; ...
                    B, sparse(Np,Np)];
            bigF = [fD; gD];
            bigu = [u; p];

            assert(~periodic);
            
            if false
                % 
                % 
                % Test with exact solution
                % 
                % 
                unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
                uI = pde.exactu(unode);
                pI = pde.exactp([node]);
                bigu = [uI(:); pI];

                %fprintf('res: \n');
                %bigA * bigu - bigF
                res_u = AD * uI(:) + BD'* pI - fD;
                RES = [[unode(:,1); unode(:,1)] [unode(:,2); unode(:,2)] AD * uI(:)  BD'* pI fD  res_u (abs(res_u) > 1e-10) * 888]
                iC = find(abs(res_u) > 1e-10)
                RES(iC, :)
                error('dbg')

                % bigFreeDof = [ufreeDof; 2*Nu+pDof];
                % bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
            end

            bigu = bigA \ bigF;
            u = bigu(1:2*Nu);
            p = bigu(2*Nu+1:end);

            residual = norm(bigF - bigA*bigu);
            info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);     
        else
            RU = rot_info.RU;
            % RU = speye(2*Nu, 2*Nu);

            %
            %  R A R' R u = R f;
            % 
            
            bigA = [RU*AD*RU', RU*BDt; ...
                    B*RU', sparse(Np,Np)];
            bigF = [RU*fD; gD];



            if ~periodic
                % set Un Dirich
                bigA(rot_info.IUn, :) = 0;
                bigA(rot_info.IUn, rot_info.IUn) = speye(length(rot_info.IUn));
                bigF(rot_info.IUn) = rot_info.Un;

                bigu = bigA \ bigF;
            else

                % 
                % reassemble
                % 
                Ndof = size(bigA, 1);
                [i,j,s] = find(bigA);
                A_period = sparse(I_period(i), I_period(j), s, Ndof, Ndof);

                F_period = zeros(Ndof, 1);
                for i = 1:Ndof
                    F_period(I_period(i)) = F_period(I_period(i)) + bigF(i);
                end


                %
                % linkage
                %
                I0 = M_period(:,1);
                I1 = M_period(:,2);
                A_period = A_period + sparse(I1, I1, ones(N_period, 1), Ndof, Ndof) + ...
                    sparse(I1, I0, -ones(N_period, 1), Ndof, Ndof);
                
                % set Un Dirich
                A_period(rot_info.IUn, :) = 0;
                A_period(rot_info.IUn, rot_info.IUn) = speye(length(rot_info.IUn));
                F_period(rot_info.IUn) = rot_info.Un;

                

                % % exact
                % if 1
                %     unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
                %     uI = pde.exactu(unode);
                %     pI = pde.exactp([node]);
                %     bigu = [RU * uI(:); pI];

                %     u_period(Idx_period(:)) = bigu(:);

                %     n1 = 2*Nu - 2*Nu_period;
                %     n2 = Np - Np_period;

                %     A11 = A_period(1:n1, 1:n1);
                %     A12 = A_period(1:n1, n2+1:n2+n1);
                %     A21 = A_period(n2+1:n2+n1, 1:n1);
                %     A22 = A_period(n2+1:n2+n1, n2+1:n2+n1);
                %     x1 = u_period(1:n1)';
                %     x2 = u_period(n2+1:n2+n1)';
                %     b1 = F_period(1:n1);
                %     b2 = F_period(n2+1:n2+n1);
                    
                %     if 1
                %         r1 = A11*x1 + A12*x2 - b1;
                %         RES = [ A11*x1 A12*x2 b1 r1 (abs(r1) > 1e-10) * 888]
                %         iC = find(abs(r1) > 1e-10)
                %         RES(iC, :)
                %     end

                %     if 0
                %         res_p = B*RU' * RU*uI(:) - g;
                %         RES = [node(:,1) node(:,2) ...
                %                B*RU' * RU*uI(:)  g  res_p (abs(res_p) > 1e-10) * 888]
                %         iC = find(abs(res_p) > 1e-10)
                %         RES(iC, :)
                %     end
                %     error('dbg')
                % end

                bigu = A_period \ F_period;
            end

            
            if  0 %true
                % 
                % 
                % Test with exact solution
                % 
                % 
                unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
                uI = pde.exactu(unode);
                pI = pde.exactp([node]);
                bigu = [RU * uI(:); pI];

                if periodic
                    A11 = A_period(1:2*Nu, 1:2*Nu);
                    A12 = A_period(1:2*Nu, 2*Nu+1:end);
                    A21 = A_period(2*Nu+1:end, 1:2*Nu);
                    A22 = A_period(2*Nu+1:end, 2*Nu+1:end);

                    b1 = F_period(1:2*Nu);
                    b2 = F_period(2*Nu+1:end);

                    u1 = bigu(1:2*Nu);
                    u2 = bigu(2*Nu+1:end);

                    if 1
                        res_u = A11*u1 + A12 * u2 - b1;
                        RES = [[unode(:,1); unode(:,1)] [unode(:,2); unode(:,2)] + slope*[unode(:,1); unode(:,1)] ...
                               A11*u1 A12*u2 b1 res_u (abs(res_u) > 1e-10) * 888]
                        iC = find(abs(res_u) > 1e-10)
                        RES(iC, :)
                    end

                    if 0
                        res_p = A21*u1 + A22 * u2 - b2;
                        RES = [node(:,1) node(:,2) + slope*node(:,1), ...
                               A21*u1 A22*u2 b2 res_p (abs(res_p) > 1e-10) * 888]
                        iC = find(abs(res_p) > 1e-10)
                        RES(iC, :)
                    end
                    
                else

                    if 1
                        res_u = RU*AD*RU' * RU*uI(:) + RU*BDt* pI - RU*fD;
                        RES = [[unode(:,1); unode(:,1)] [unode(:,2); unode(:,2)] + slope*[unode(:,1); unode(:,1)] ...
                               RU*AD*RU' * RU*uI(:)  RU*BDt* pI  RU*fD  res_u (abs(res_u) > 1e-10) * 888]
                        iC = find(abs(res_u) > 1e-10)
                        RES(iC, :)
                    end

                    if 0
                        res_p = B*RU' * RU*uI(:) - g;
                        RES = [node(:,1) node(:,2) ...
                               B*RU' * RU*uI(:)  g  res_p (abs(res_p) > 1e-10) * 888]
                        iC = find(abs(res_p) > 1e-10)
                        RES(iC, :)
                    end
                end
                
                
                error('dbg')

                % bigFreeDof = [ufreeDof; 2*Nu+pDof];
                % bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
            end
            
            u_period = bigu(1:2*Nu);
            p = bigu(2*Nu+1:end);

            u = RU' * u_period;

            %error('period check');
        end
        

      case 'mg'
        %         option.tol = Np^(-2);        
        option.solver  = 'WCYCLE';
        [u(ufreeDof),p,info] = mgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                        u(ufreeDof),p,elem,ufreeDof,option);         
      case 'asmg'
        [u(ufreeDof),p,info] = asmgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                          u,p,node,elem,bdFlag,ufreeDof,option); 
    end

    %% Post-process
    if length(pDof) ~= Np % p is unique up to a constant
                          % impose the condition int(p)=0
        c = sum(mean(p(elem),2).*area)/sum(area);
        p = p - c;
    end

    %% Output
    soln = struct('u',u,'p',p);
    eqn = struct('A',AD,'Bt',BDt,'Lap',A,'f',f,'g',g,...
                 'edge',edge,'ufreeDof',ufreeDof,'pDof',pDof, ...
                 'A_period', A_period, 'F_period', F_period, 'sol_period', bigu);
    info.assembleTime = assembleTime;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subfunctions getbdStokesP2P1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,BDt,fD,gD,u,p,ufreeDof,fixedDof,pDof,rot_info] = getbdStokesP2P1
    %% Boundary condition of Stokes equation: P2-P0 elements

    %% Initial set up
    % f = [f1; f2];    % set in Neumann boundary condition
    % g = zeros(Np,1);
        u = zeros(2*Nu,1);    
        p = zeros(Np,1);
        ufreeDof = (1:Nu)';
        pDof = (1:Np)';
        rot_info = [];
        
        if ~exist('bdFlag','var'), bdFlag = []; end
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'beta')
            if isfield(pde,'g_R')
                pde.beta = pde.g_R;
            else
                pde.beta = [];
            end
        end
        if ~isfield(pde,'g_R')
            if isfield(pde,'g_RN')
                pde.g_R = pde.g_RN;
            else
                pde.g_R = [];
            end
        end
        if isfield(pde,'g_Dn')
            gDNormal = pde.g_Dn;
        else
            gDNormal = pde.g_D;
        end
               
        %% Part 1: Find Dirichlet dof and modify the matrix
        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedDof = false(Nu,1);     
        if ~isempty(bdFlag)       % case: bdFlag is not empty 
            elem2edge = elem2dof(:,4:6)-N;
            isDirichlet(elem2edge(bdFlag(:)==1)) = true;
            isFixedDof(edge(isDirichlet,:)) = true;   % nodes of all D-edges
            isFixedDof(N + find(isDirichlet')) = true;% dof on D-edges
            fixedDof = find(isFixedDof);
            ufreeDof = find(~isFixedDof);            
        end
        if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.beta)
            fixedDof = bdDof; 
            isFixedDof(fixedDof) = true;
            ufreeDof = find(~isFixedDof);    
        end

        % if isempty(fixedDof) % pure Neumann boundary condition
        %                      % pde.g_N could be empty which is homogenous Neumann boundary condition
        %     fixedDof = 1;
        %     ufreeDof = (2:Nu)';    % eliminate the kernel by enforcing u(1) = 0;
        % end

        % Modify the matrix
        % Build Dirichlet boundary condition into the matrix AD by enforcing
        % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
        % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.

        % % disabled
        % bdidx = zeros(2*Nu,1); 
        % bdidx([fixedDof; Nu+fixedDof]) = 1;
        % Tbd = spdiags(bdidx,0,2*Nu,2*Nu);
        % T = spdiags(1-bdidx,0,2*Nu,2*Nu);
        % AD = T*A*T + Tbd;
        % BD = B*T;

        

        %% Part 2: Find boundary edges and modify the right hand side f and g
        % Find boundary edges: Neumann and Robin
        Neumann = []; Robin = []; %#ok<*NASGU>
        if ~isempty(bdFlag)
            isNeumann(elem2edge(bdFlag(:)==2)) = true;
            NeumannIdx = find(isNeumann);        
            Neumann   = edge(isNeumann,:);
            nNeumann = size(Neumann, 1);
            %Neumann
            %NeumannIdx

            isRobin(elem2edge(bdFlag(:)==3)) = true;
            Robin     = edge(isRobin,:);
            RobinIdx = find(isRobin);        
            Robin   = edge(isRobin,:);
            nRobin = size(Robin, 1);
            % Robin
            % RobinIdx
            % error('dbg')

            if use_slip
                % build rotation

                % rotation per elem
                ve = node(Robin(:,1),:) - node(Robin(:,2),:);
                ne = [-ve(:,2), ve(:,1)] ./ sqrt(ve(:,1).^2 + ve(:,2).^2);
                iflip = find(ne(:,2) > 0);
                ne(iflip,:) = -ne(iflip,:);
                te = [-ne(:,2), ne(:,1)];
                
                % rotation per node
                if 0
                    % for general case
                    ne_sum = zeros(N, 2);
                    ne_sum(:, 1) = 1;

                    ne_sum(Robin(:)) = 0;
                    ne_sum(Robin(:, 1), :) = ne(:, :);
                    ne_sum(Robin(:, 2), :) = ne_sum(Robin(:, 2), :) + ne(:, :);
                    ne_sum(unique(Robin(:)), :);
                    IR = unique(Robin(:));
                    
                    for i = 1:length(IR)
                        % normalize
                        ii = IR(i);
                        norm(ne_sum(ii, :));
                        ne_sum(ii, :) = ne_sum(ii, :) ./ norm(ne_sum(ii, :));
                    end
                    ne_sum(IR, :);

                    RU = speye(2*Nu, 2*Nu);
                    % RU(IR, IR) = 0;
                    % RU(Nu+IR, Nu+IR) = 0;

                    RU(1:N, 1:N)                 = diag(ne_sum(:, 1));
                    RU(1:N, Nu+1:Nu+N)          = diag(ne_sum(:, 2));
                    RU(Nu+1:Nu+N, 1:N)           = diag(-ne_sum(:, 2));
                    RU(Nu+1:Nu+N, Nu+1:Nu+N)     = diag(ne_sum(:, 1));

                    % RU(IR, IR)                   = ne_sum(:, 1);
                    % RU(IR, Nu+IR)                = ne_sum(:, 2);
                    % RU(Nu+IR, IR)                = ne_sum(:, 2);
                    % RU(Nu+IR, Nu+IR)             = -ne_sum(:, 1);

                    RU(N+RobinIdx, N+RobinIdx)       = diag(ne(:, 1));
                    RU(N+RobinIdx, Nu+N+RobinIdx)    = diag(ne(:, 2));
                    RU(Nu+N+RobinIdx, N+RobinIdx)    = diag(-ne(:, 2));
                    RU(Nu+N+RobinIdx, Nu+N+RobinIdx) = diag(ne(:, 1));

                else
                    % square domain only
                    ne_sum = zeros(nRobin + 1, 2);
                    ne_sum(1:nRobin,   :) = ne(:, :);
                    ne_sum(2:nRobin+1, :) = ne_sum(2:nRobin+1, :) + ne(:, :);

                    if periodic
                        ne_sum(nRobin+1, :) = ne_sum(nRobin+1, :) + ne(1, :);
                        ne_sum(1, :) = ne_sum(1, :) + ne(nRobin, :);
                    end
                    
                    for i = 1:nRobin+1
                        % normalize
                        ne_sum(i, :) = ne_sum(i, :) ./ norm(ne_sum(i, :));
                    end


                    RU = speye(2*Nu, 2*Nu);

                    IR = [Robin(:, 1); Robin(end, 2)];

                    RU(IR,    IR)    = diag(ne_sum(:,  1));
                    RU(IR,    Nu+IR) = diag(ne_sum(:,  2));
                    RU(Nu+IR, IR)    = diag(-ne_sum(:, 2));
                    RU(Nu+IR, Nu+IR) = diag(ne_sum(:,  1));
                    
                    RU(N+RobinIdx, N+RobinIdx)       = diag(ne(:, 1));
                    RU(N+RobinIdx, Nu+N+RobinIdx)    = diag(ne(:, 2));
                    RU(Nu+N+RobinIdx, N+RobinIdx)    = diag(-ne(:, 2));
                    RU(Nu+N+RobinIdx, Nu+N+RobinIdx) = diag(ne(:, 1));
                    
                    %error('dbg ne_sum')
                end

                

                % 
                % Note: !!! 
                %  
                % Rotate matrix act on periodic dofs, 
                %    but do not apply dirich B.C. on the duplicate dofs, 
                %    instead maintain linkage between peridoic dofs.
                % 
                % 
                if ~periodic
                    IUn = [IR; N+RobinIdx';];
                else
                    IUn = [IR(1:nRobin); N+RobinIdx';];
                end

                rot_info.RU = RU;
                rot_info.IUn = IUn;

                if isnumeric(gDNormal)
                    if ~periodic
                        rot_info.Un = gDNormal([1:nRobin, 1, nRobin+1:2*nRobin]);
                    else
                        rot_info.Un = gDNormal([1:nRobin, nRobin+1:2*nRobin]);
                    end
                else
                        rot_info.Un = gDNormal([node(IR,:); 
                                                (node(edge(RobinIdx',1),:)+node(edge(RobinIdx',2),:))/2]);
                end
                
                clear ne_sum ne ve te;
                %error('dbg')
            end
        end

        % if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.beta))
        %     % no bdFlag, only pde.g_N or pde.beta is given in the input
        %     Neumann = edge(bdDof>N,:);
        %     if ~isempty(pde.beta)
        %         Robin = Neumann;
        %     end
        % end

        % Neumann boundary condition
        if ~isempty(Neumann) 
            [lambda,w] = quadpts1(quadorder1d);
            nQuad = size(lambda,1);

            % quadratic bases (1---3---2)
            bdphi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
            bdphi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
            bdphi(:,3) = 4*lambda(:,1).*lambda(:,2);

            % length of edge
            ve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
            edgeLength = sqrt(sum(ve.^2,2));

            % update RHS
            gex = zeros(size(Neumann,1),2);   % x-component
            gey = zeros(size(Neumann,1),2);   % y-component
            for pp = 1:nQuad
                pxy = lambda(pp,1)*node(Neumann(:,1),:)+lambda(pp,2)*node(Neumann(:,2),:);
                
                if isnumeric(pde.g_N)
                    assert(size(pde.g_N, 1) == 2 * nNeumann);

                    gp1 = pde.g_N([1:nNeumann],:);
                    gp2 = pde.g_N([2:nNeumann, 1],:);
                    gp3 = pde.g_N([nNeumann+1:2*nNeumann],:);
                    gp = gp1 .*  bdphi(pp,1) ...
                         + gp2 .*  bdphi(pp,2) ...
                         + gp3 .*  bdphi(pp,3);

                    if size(pde.g_N, 2) == 2
                        ;
                    elseif size(pde.g_N, 2) == 4
                        % g_N = f1 * f2
                        gp(:, 1) = gp(:, 1).*gp(:, 2);
                        gp(:, 2) = gp(:, 3).*gp(:, 4);
                        gp = gp(:, [1 2]);
                    else
                        error('g_N dim not right')
                    end
                    
                else
                    gp = pde.g_N(pxy);
                end
                
                gex(:,1) = gex(:,1) + w(pp)*edgeLength.*gp(:,1)*bdphi(pp,1);
                gex(:,2) = gex(:,2) + w(pp)*edgeLength.*gp(:,1)*bdphi(pp,2);
                gey(:,1) = gey(:,1) + w(pp)*edgeLength.*gp(:,2)*bdphi(pp,1);
                gey(:,2) = gey(:,2) + w(pp)*edgeLength.*gp(:,2)*bdphi(pp,2);
                f1(N+NeumannIdx) = f1(N+NeumannIdx) + w(pp)*edgeLength.*gp(:,1)*bdphi(pp,3); % interior bubble
                f2(N+NeumannIdx) = f2(N+NeumannIdx) + w(pp)*edgeLength.*gp(:,2)*bdphi(pp,3); % interior bubble
            end
            f1(1:N) = f1(1:N) + accumarray(Neumann(:), gex(:),[N,1]);
            f2(1:N) = f2(1:N) + accumarray(Neumann(:), gey(:),[N,1]);
            % unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
            % [[unode(:,1); unode(:,1)] [unode(:,2); unode(:,2)] [f1; f2]]
        end
        
        % The case non-empty Neumann but g_N=[] corresponds to the zero flux
        % boundary condition on Neumann edges and no modification is needed.





        % Robin boundary condition
        if ~isempty(Robin) 
            [lambda,w] = quadpts1(quadorder1d);
            nQuad = size(lambda,1);

            % quadratic bases (1---3---2)
            bdphi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
            bdphi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
            bdphi(:,3) = 4*lambda(:,1).*lambda(:,2);

            % length of edge
            ve = node(Robin(:,1),:) - node(Robin(:,2),:);
            edgeLength = sqrt(sum(ve.^2,2));

            % update RHS
            gex = zeros(size(Robin,1),2);   % x-component
            gey = zeros(size(Robin,1),2);   % y-component
            ss = zeros(size(Robin,1),3,3);

            xqt = [];
            mqt = [];
            
            for pp = 1:nQuad
                pxy = lambda(pp,1)*node(Robin(:,1),:)+lambda(pp,2)*node(Robin(:,2),:);

                if isnumeric(pde.beta)
                    assert(size(pde.beta, 1) == 2 * nRobin);
                    
                    gr1 = pde.beta([1:nRobin], :);
                    gr2 = pde.beta([2:nRobin, 1], :);
                    gr3 = pde.beta([nRobin+1:2*nRobin], :);
                    gr = gr1 .*  bdphi(pp,1) ...
                         + gr2 .*  bdphi(pp,2) ...
                         + gr3 .*  bdphi(pp,3);
                    % xqt = [xqt pxy(:,1)];
                    % mqt = [mqt gr];

                    % gr relates to basis phi_i
                    % if use_linear
                    %     gr = [pde.beta([1:nRobin]), pde.beta([2:nRobin, 1]), .5*(pde.beta([1:nRobin]) + pde.beta([2:nRobin, 1]))];
                    % else

                    % gr = [pde.beta([1:nRobin]), pde.beta([2:nRobin, 1]), pde.beta([nRobin+1:2*nRobin])];

                    %end

                    gp1 = pde.g_R([1:nRobin], :);
                    gp2 = pde.g_R([2:nRobin, 1], :);
                    gp3 = pde.g_R([nRobin+1:2*nRobin], :);
                    
                    gp = gp1 .*  bdphi(pp,1) ...
                         + gp2 .*  bdphi(pp,2) ...
                         + gp3 .*  bdphi(pp,3);

                    if size(pde.g_R, 2) == 2
                        ;
                    elseif size(pde.g_R, 2) == 4
                        % g_R = f1 * f2
                        gp(:, 1) = gp(:, 1).*gp(:, 2);
                        gp(:, 2) = gp(:, 3).*gp(:, 4);
                        gp = gp(:, [1 2]);
                    else
                        error('g_R dim not right')
                    end
                    
                else
                    gp = pde.g_R(pxy);
                    gr = pde.beta(pxy);
                end

                % 1st component
                gex(:,1) = gex(:,1) + w(pp)*edgeLength.*gp(:,1)*bdphi(pp,1); % to Robin(:,1)
                gex(:,2) = gex(:,2) + w(pp)*edgeLength.*gp(:,1)*bdphi(pp,2); % to Robin(:,2)
                f1(N+RobinIdx) = f1(N+RobinIdx) + w(pp)*edgeLength.*gp(:,1)*bdphi(pp,3); % interior bubble, 
                                                                                         % to N + RobinIdx

                % 2nd component
                gey(:,1) = gey(:,1) + w(pp)*edgeLength.*gp(:,2)*bdphi(pp,1);
                gey(:,2) = gey(:,2) + w(pp)*edgeLength.*gp(:,2)*bdphi(pp,2);
                f2(N+RobinIdx) = f2(N+RobinIdx) + w(pp)*edgeLength.*gp(:,2)*bdphi(pp,3); % interior bubble

                for iR = 1:3
                    for jR = 1:3   
                        ss(:,iR,jR) = ss(:,iR,jR) + ...
                            w(pp)*edgeLength.*gr(:)*bdphi(pp,iR)*bdphi(pp,jR);
                        %weightgR(pp)*gRp*phigR(pp,iR).*phigR(pp,jR);
                    end
                end
                
            end
            f1(1:N) = f1(1:N) + accumarray(Robin(:), gex(:),[N,1]);
            f2(1:N) = f2(1:N) + accumarray(Robin(:), gey(:),[N,1]);
            % unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
            % [[unode(:,1); unode(:,1)] [unode(:,2); unode(:,2)] [f1; f2]]

            R1 = Robin(:,1);
            R2 = Robin(:,2);
            R3 = N+RobinIdx';
            jj = [R1; R1; R1; R2; R2; R2; R3; R3; R3];
            ii = [R1; R2; R3; R1; R2; R3; R1; R2; R3];
            AR = sparse(ii, jj, ss(:), Nu, Nu);

            A = A + blkdiag(AR,AR);

            if dbg_on
                data = [xqt(:) mqt(:)];
                data = sortrows(data);
                plot(data(:,1), data(:,2), '-o');
                error('stop')
            end
        end






        
        % Dirichlet boundary conditions
        fD = [f1; f2];
        %norm(fD)
        if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
            u1 = zeros(Nu,1);
            u2 = zeros(Nu,1);
            idx = (fixedDof > N);              % index of edge dof
            uD = pde.g_D(node(fixedDof(~idx),:));  % bd value at vertex dofs    
            u1(fixedDof(~idx)) = uD(:,1);
            u2(fixedDof(~idx)) = uD(:,2);

            bdEdgeIdx = fixedDof(idx)-N;
            bdEdgeMid = (node(edge(bdEdgeIdx,1),:)+node(edge(bdEdgeIdx,2),:))/2;
            uD = pde.g_D(bdEdgeMid);         % bd values at middle points of edges
            u1(fixedDof(idx)) = uD(:,1);
            u2(fixedDof(idx)) = uD(:,2);

            % u = [u1; u2]; % Dirichlet bd condition is built into u
            % f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
            % g = g - B*u;  % the right hand side
            % g = g - mean(g);
            
            fD(fixedDof)    = u1(fixedDof);
            fD(fixedDof+Nu) = u2(fixedDof);
        end
        %norm(fD)
        % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
        % boundary condition and no modification is needed.
        
        % modfiy pressure dof for pure Dirichlet
        if isempty(Neumann)
            pDof = (1:Np-1)';
        end
        
        ufreeDof = [ufreeDof; Nu+ufreeDof];


        % Dirichlet
        AD = A;
        BDt = B';

        if option.verb > 1
            fprintf('fixed: \n');
            fixedDof
        end
        
        bdidx = zeros(2*Nu,1); 
        bdidx([fixedDof; Nu+fixedDof]) = 1;
        Tbd = spdiags(bdidx,0,2*Nu,2*Nu);
        AD([fixedDof; Nu+fixedDof], :) = 0;
        AD = AD + Tbd;
        BDt([fixedDof; Nu+fixedDof], :) = 0;

        gD = g;
        
    end % end of function getbdStokesP2P1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
