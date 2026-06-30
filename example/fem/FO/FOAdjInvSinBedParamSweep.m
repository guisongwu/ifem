%% FOADJINVSINBEDPARAMSWEEP Parameter sweep for FO SinBed inversion.
%
% Compare lambda, alpha, and eps_reg for FOAdjInvSinBed.  The main metrics
% are inverse iterations, final objective, data residual, and beta Linf
% error.  Edit the lists below to widen or narrow the scan.

close all;
clear variables;

alphaList = [0,1e-13,1e-12,1e-11];
lambdaList = [1e-8,1e-7];
epsRegList = [1e-4,1e-3];

maxInverseIt = 50;
useLineSearch = false;

nCase = numel(alphaList)*numel(lambdaList)*numel(epsRegList);
result(nCase,1) = struct(...
    'alpha',NaN,...
    'lambda',NaN,...
    'eps_reg',NaN,...
    'iterations',NaN,...
    'objective',NaN,...
    'dataResidual',NaN,...
    'betaLinfAbs',NaN,...
    'betaLinfRel',NaN,...
    'gradientNorm',NaN,...
    'forwardSolves',NaN,...
    'failed',false,...
    'message','');

caseId = 0;
fprintf(['FO SinBed parameter sweep: %d cases, ',...
    'maxInverseIt = %d, lineSearch = %d\n'],...
    nCase,maxInverseIt,useLineSearch);
fprintf('%3s %10s %10s %10s %4s %12s %12s %12s %12s %12s %6s\n',...
    'id','alpha','lambda','eps_reg','it','objective',...
    'dataRes','betaLinf','betaRel','|grad|','fSolve');

for ie = 1:numel(epsRegList)
    for il = 1:numel(lambdaList)
        for ia = 1:numel(alphaList)
            caseId = caseId+1;
            foInvConfig = struct(...
                'alpha',alphaList(ia),...
                'lambda',lambdaList(il),...
                'eps_reg',epsRegList(ie),...
                'maxInverseIt',maxInverseIt,...
                'useLineSearch',useLineSearch,...
                'figureVisible','off');

            result(caseId).alpha = foInvConfig.alpha;
            result(caseId).lambda = foInvConfig.lambda;
            result(caseId).eps_reg = foInvConfig.eps_reg;

            try
                evalc('run(''FOAdjInvSinBed.m'')');
                result(caseId).iterations = numel(history.objective);
                result(caseId).objective = history.objective(end);
                result(caseId).dataResidual = history.dataResidual(end);
                result(caseId).betaLinfAbs = betaErrorLinf;
                result(caseId).betaLinfRel = betaErrorRelativeLinf;
                result(caseId).gradientNorm = history.gradientNorm(end);
                result(caseId).forwardSolves = optimizationForwardSolves;
            catch exception
                result(caseId).failed = true;
                result(caseId).message = exception.message;
            end

            printresultrow(caseId,result(caseId));
        end
    end
end

result = result(1:caseId);
failed = arrayfun(@(item) item.failed,result);
valid = find(~failed);
if isempty(valid)
    fprintf('\nAll parameter-sweep cases failed.\n');
else
    betaLinfAbs = arrayfun(@(item) item.betaLinfAbs,result(valid));
    [~,order] = sort(betaLinfAbs);
    sorted = valid(order);

    fprintf('\nBest cases sorted by beta Linf error\n');
    fprintf('%3s %10s %10s %10s %4s %12s %12s %12s %12s %12s %6s\n',...
        'id','alpha','lambda','eps_reg','it','objective',...
        'dataRes','betaLinf','betaRel','|grad|','fSolve');
    for k = 1:min(10,numel(sorted))
        printresultrow(sorted(k),result(sorted(k)));
    end

    best = result(sorted(1));
    fprintf(['\nBest by beta Linf: alpha = %.3e, lambda = %.3e, ',...
        'eps_reg = %.3e, iterations = %d, objective = %.6e, ',...
        'betaLinfAbs = %.6e, betaLinfRel = %.6e\n'],...
        best.alpha,best.lambda,best.eps_reg,best.iterations,...
        best.objective,best.betaLinfAbs,best.betaLinfRel);
end

function printresultrow(caseId,item)
    if item.failed
        fprintf('%3d %10.3e %10.3e %10.3e failed: %s\n',...
            caseId,item.alpha,item.lambda,item.eps_reg,item.message);
    else
        fprintf('%3d %10.3e %10.3e %10.3e %4d %12.4e %12.4e %12.4e %12.4e %12.4e %6d\n',...
            caseId,item.alpha,item.lambda,item.eps_reg,...
            item.iterations,item.objective,item.dataResidual,...
            item.betaLinfAbs,item.betaLinfRel,item.gradientNorm,...
            item.forwardSolves);
    end
end
