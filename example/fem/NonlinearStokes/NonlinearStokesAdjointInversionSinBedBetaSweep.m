%% NONLINEARSTOKESADJOINTINVERSIONSINBEDBETASWEEP Geometry sensitivity sweep.
%
% Run several sinusoidal-bed beta inversions with different (L,H).  The
% mesh is resized with the physical domain so that dx and dz stay close to
% targetMeshSize.  Boundary flags are set on the reference strip
% [0,L]x[0,1], where the top is always y=1 before the B-type mapping.

close all;
clear variables;
set(groot,'DefaultFigureVisible','off');

targetMeshSize = 0.1;
caseList = [
    0.5, 0.5
    1.0, 0.5
    2.0, 0.5
    1.0, 1.0
    2.0, 1.0
    4.0, 1.0
    ];

nCase = size(caseList,1);
result = table('Size',[nCase,13],...
    'VariableTypes',{'double','double','double','double','double',...
                     'double','logical','double','double','double',...
                     'double','double','double'},...
    'VariableNames',{'L','H','Nx','Nz','Nm','nObs','converged',...
                     'iterations','objective','dataResidual',...
                     'betaL2Relative','betaLinfRelative','gradientNorm'});

for iCase = 1:nCase
    L = caseList(iCase,1);
    H = caseList(iCase,2);
    Nx = max(4,round(L/targetMeshSize));
    Nz = max(3,round(H/targetMeshSize));
    caseResult = runinversioncase(L,H,Nx,Nz);
    result.L(iCase) = L;
    result.H(iCase) = H;
    result.Nx(iCase) = Nx;
    result.Nz(iCase) = Nz;
    result.Nm(iCase) = caseResult.Nm;
    result.nObs(iCase) = caseResult.nObs;
    result.converged(iCase) = caseResult.converged;
    result.iterations(iCase) = caseResult.iterations;
    result.objective(iCase) = caseResult.objective;
    result.dataResidual(iCase) = caseResult.dataResidual;
    result.betaL2Relative(iCase) = caseResult.betaL2Relative;
    result.betaLinfRelative(iCase) = caseResult.betaLinfRelative;
    result.gradientNorm(iCase) = caseResult.gradientNorm;

    fprintf(['L=%4.1f H=%4.1f Nx=%2d Nz=%2d Nm=%2d nObs=%2d ',...
        'it=%2d obj=%.04e betaL2=%.04e betaLinf=%.04e\n'],...
        L,H,Nx,Nz,caseResult.Nm,caseResult.nObs,...
        caseResult.iterations,caseResult.objective,...
        caseResult.betaL2Relative,caseResult.betaLinfRelative);
end

disp(result);

function result = runinversioncase(L,H,Nx,Nz)
slope = tan(0.5*pi/180);
bedAmplitude = 0.1*H;

[node,elem] = rectanglemesh(L,1,Nx,Nz);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node = maptoexperimentb(node,L,H,bedAmplitude,slope);

[~,edge] = dofP2(elem);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
surfaceLevel = -slope*uNode(:,1);
tolGeometry = 100*eps(max(1,max(abs(node(:)))));
topDof = find(abs(uNode(:,2)-surfaceLevel)<tolGeometry ...
            & uNode(:,1)<L-tolGeometry);
[~,order] = sort(uNode(topDof,1));
topDof = topDof(order);
xObs = uNode(topDof,1);
topWeight = boundaryweights(xObs,L,slope);

pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-11;
option.residual_tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
option.assemble_tangent = true;

Nm = Nx;
xBeta = (0:Nm-1)'*L/Nm;
betaTrue = 2*(1+0.25*cos(2*pi*xBeta/L));
betaInitial = betaTrue.*(1+0.2*sin(2*pi*xBeta/L)+0.05);
qTrue = log(betaTrue);
q = log(betaInitial);

[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(trueInfo.converged,'Truth solve failed for L=%g, H=%g.',L,H);
dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

maxInverseIt = 20;
lambda = 1e-7;
pcgTolerance = 1e-8;
pcgMaxIt = 50;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;
uWarm = [];
converged = false;

for k = 1:maxInverseIt
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    assert(forwardInfo.converged,...
        'Forward solve failed for L=%g, H=%g at iteration %d.',L,H,k);
    uWarm = u;
    residual = u(topDof)-dataObs;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    dataResidual = sqrt((topWeight'*(residual.^2))/dataNormSquared);
    G = assembleparameterderivative(eqn,q,xBeta,L,Nm);

    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint;

    if norm(gradient) <= gradientTolerance
        converged = true;
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,lambda);
    [step,~,~,~] = pcg(hessian,-gradient,pcgTolerance,pcgMaxIt);
    if norm(step) <= stepTolerance*max(1,norm(q))
        converged = true;
        break
    end

    accepted = false;
    stepLength = 1;
    for lineSearchIt = 1:10
        qTrial = q+stepLength*step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,L);
        if trialInfo.converged
            trialResidual = uTrial(topDof)-dataObs;
            trialObjective = 0.5*(topWeight'*(trialResidual.^2)) ...
                /dataNormSquared;
            if trialObjective < objective
                q = qTrial;
                uWarm = uTrial;
                lambda = max(lambda/3,1e-12);
                accepted = true;
                break
            end
        end
        stepLength = stepLength/2;
    end
    if ~accepted
        lambda = 10*lambda;
    end
end

betaRecovered = exp(q);
result.Nm = Nm;
result.nObs = numel(topDof);
result.converged = converged;
result.iterations = k;
result.objective = objective;
result.dataResidual = dataResidual;
result.betaL2Relative = norm(betaRecovered-betaTrue)/norm(betaTrue);
result.betaLinfRelative = norm(betaRecovered-betaTrue,inf)/norm(betaTrue,inf);
result.gradientNorm = norm(gradient);
end

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,lambda)
incrementalState = eqn.tangent\(-G*direction);
incrementalObservation = zeros(size(eqn.tangent,1),1);
incrementalObservation(topDof) = ...
    topWeight.*incrementalState(topDof)/dataNormSquared;
incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
product = G'*incrementalAdjoint+lambda*direction;
end

function G = assembleparameterderivative(eqn,q,xBeta,L,Nm)
G = zeros(size(eqn.tangent,1),Nm);
beta = exp(q(:));
for j = 1:Nm
    direction = zeros(Nm,1);
    direction(j) = 1;
    deltaBeta = beta.*direction;
    directionFunction = @(pt) periodicP1(...
        pt(:,1),xBeta,deltaBeta,L);
    G(:,j) = eqn.applyBetaDerivative(directionFunction);
end
end

function [u,eqn,info] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L)
beta = exp(q(:));
pde.beta = @(pt) periodicP1(pt(:,1),xBeta,beta,L);
if isempty(u0)
    if isfield(option,'u0')
        option = rmfield(option,'u0');
    end
else
    option.u0 = u0;
end
[soln,eqn,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option);
u = soln.u;
end

function value = periodicP1(x,xNode,nodalValue,L)
xWrapped = mod(x,L);
value = interp1([xNode;L],[nodalValue;nodalValue(1)],...
    xWrapped,'linear');
end

function weight = boundaryweights(xObs,L,slope)
nObs = numel(xObs);
assert(nObs > 0,'No top-boundary observation dofs were found.');
weight = sqrt(1+slope^2)*(L/nObs)*ones(nObs,1);
end

function node = maptoexperimentb(node,L,H,bedAmplitude,slope)
x = node(:,1);
r = node(:,2);
surface = -slope*x;
bed = surface-H+bedAmplitude*sin(2*pi*x/L);
node(:,2) = bed+r.*(surface-bed);
end

function [node,elem] = rectanglemesh(L,H,Nx,Nz)
x = linspace(0,L,Nx+1);
z = linspace(0,H,Nz+1);
[X,Z] = meshgrid(x,z);
node = [X(:),Z(:)];

cellId = reshape(1:(Nx+1)*(Nz+1),Nz+1,Nx+1);
elem = zeros(2*Nx*Nz,3);
cursor = 0;
for ix = 1:Nx
    for iz = 1:Nz
        v1 = cellId(iz,ix);
        v2 = cellId(iz,ix+1);
        v3 = cellId(iz+1,ix);
        v4 = cellId(iz+1,ix+1);
        cursor = cursor+1;
        elem(cursor,:) = [v1,v2,v4];
        cursor = cursor+1;
        elem(cursor,:) = [v1,v4,v3];
    end
end
end
