%% Test Poisson

close all; 
clear variables;

maxIt = 1; 
errL2 = zeros(maxIt,1); 

pde = test_P1_data;
%pde = test_P2_data;

h = 0.1;

for k = 1:maxIt

    [node,elem] = squaremesh([0 1 0 1], h);
    bdFlag = setboundary(node,elem,'Dirichlet');

    [soln,eqn,info] = Poisson(node,elem,bdFlag,pde);
    uh = soln.u;

    figure(1);
    showresult(node,elem,uh);

    errL2(k) = getL2error(node,elem,pde.exactu,uh);
    fprintf('pass: %d, h: %f, err: %e\n', k, h, errL2(k));

    h = h / 2;
end


if maxIt > 1
    fprintf('\nConvergence rate\n');
    for k = 1:maxIt
        if k == 1
            fprintf('err: %e -\n', errL2(k));
        else
            fprintf('err: %e %.2f\n', errL2(k), -log(errL2(k)/errL2(k-1)) / log(2) );
        end
    end
end
