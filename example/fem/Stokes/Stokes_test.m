%% SQUARESTOKE Stokes equations on the unit square
%
%   SQUARESTOKE computes P2-P1 approximations of the Stokes equations in
%   the unit square on a sequence of meshes obtained by uniform refinement.
%   It plots the approximation error (pressue in L2 norm and velocity in H1
%   norm) vs the number of dof. Other types of FEM approximation can be
%   computed similarly.
% 
% See also StokesP2P1, collidingflow, squarePoisson 
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
clear variables;


%% Set up
maxIt = 1;
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
erru = zeros(maxIt,1); 
errp = zeros(maxIt,1);




%% PDE and options
% pde = Stokesdata2;
% pde = Stokesdata3;
pde = test_Stokes_P1_data;

h = 0.05;
option.solver = 'direct';

%% Finite Element Method        
for k = 1:maxIt

    [node,elem] = squaremesh([0 1 0 1], h);
    bdFlag = setboundary(node,elem,'Dirichlet');
    [soln,eqn] = StokesP2P1(node,elem,bdFlag,pde,option);

    nnode = size(node, 1);
    ndofP2 = size(node, 1) + size(eqn.edge, 1); % P2
    
    uh = soln.u;
    ph = soln.p;
    N(k) = length(uh)+length(ph);
    h(k) = 1./(sqrt(size(node,1))-1);

    if k == 1
        figure(1);
        showresult(node,elem,uh(ndofP2+1:2*ndofP2));
    end

    % compute error
    uI = pde.exactu([node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2]);
    % erru(k) = sqrt((uh-uI(:))'*eqn.A*(uh-uI(:)));
    
    erru(k) = getL2error(node,elem,pde.exactu,reshape(uh, ndofP2, 2));
    errp(k) = getL2error(node,elem,pde.exactp,ph);
    fprintf('pass: %d, h: %f, err: %e %e\n', k, h, erru(k), errp(k));

    h = h / 2;
end



%% Plot convergence rates
if maxIt > 1
    fprintf('\nConvergence rate\n');

    for k = 1:maxIt
        if k == 1
            fprintf('err_u: %e -\n', erru(k));
        else
            fprintf('err_u: %e %.2f\n', erru(k), -log(erru(k)/erru(k-1)) / log(2) );
        end
    end
    fprintf('\n');
    for k = 1:maxIt
        if k == 1
            fprintf('err_p: %e -\n', errp(k));
        else
            fprintf('err_p: %e %.2f\n', errp(k), -log(errp(k)/errp(k-1)) / log(2) );
        end
    end
    fprintf('\n');
end


