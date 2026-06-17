%% Stokes inversion for m

fprintf('==============Stokes inversion robin beta===============\n');
global slope h dbg_case dbg_on;                           % slab slope

dbg_on = false;
dbg_case = 0;

%% Control parameters
slope = 0.1;
h = 0.1;
scheme = 6;                             
max_iteration = 5;

% plot dim
plot_m = 2;
plot_n = 3;
option.verb = 0;
option.solver = 'direct';
option.quadorder = 4;
option.use_newton = false;
option.use_slip = true;
option.periodic = true;


%% Constants
four = mp('4');
three = mp('3');
two = mp('2');
one = mp('1');


%% Setup domain
figure(1);
[node,elem] = squaremesh([0 1 0 0.5], h);
bdFlag = setboundary(node,elem, 'Neumann','y==0.5', 'Robin', 'y==0');
nx = length([0:h:1]);
ny = length([0:h:.5]);
n1 = nx-1;

X = reshape(node(:,1), ny, nx);
Y = reshape(node(:,2), ny, nx);
assert(norm(X(1,:)' - [0:h:1]') < 1e-12);
assert(norm(Y(:,1) - [0:h:.5]') < 1e-12);

N = size(node, 1);
[elem2dof,edge,bdDof] = dofP2(elem);
NE = size(edge,1);
Nu = N + NE;    % velocity u using P2
Np = N;         % pressure p using P1 
Nm = n1;        % parameter m with n1 dof
EI = eye(Nm);

if slope ~= 0
    %
    % Slab test remap coordinate
    % Note: this must be done after bdFlag is set
    % 
    fprintf(2, 'Slab test %f\n', slope);
    node(:,2) = node(:,2) - slope * node(:,1);
end

%% Index
IUxNode = [1:N];                 % extract u node value from uv
IUyNode = [Nu+1:Nu+N];           % extract v node value from uv

unode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
unode(:,2) = unode(:,2) + slope * unode(:,1);

IUxBot = sort(find(abs(unode(:,2)) < 1e-8 & unode(:,1) < 1-h/4));
IUyBot = IUxBot + Nu;

IUxTop = sort(find(abs(unode(:,2) - .5) < 1e-8 & unode(:,1) < 1-h/4));
IUyTop = IUxTop + Nu;

Isft = reshape([1:n1; n1+1:2*n1], 2*n1, 1); 
sft = @(dat) dat(Isft);
sft_ext = @(dat) sft(extend_mid(dat));
 

%% Setup pde and solve 
pde = stokes_data_grav_period;
warning('\nChange function to data !!!\n');

xbot = [0:h:1-h h/2:h:1]';              % periodic without last one
                                        % first node, then edge center
ybot = 0 - slope * xbot;
pt_bot = [xbot, ybot];

% top
xtop = xbot;
ytop = 0.5 - slope * xbot;
pt_top = [xbot, ytop];

% bot
pde.g_N = pde.g_N(pt_top);
pde.g_RN = pde.g_RN(pt_bot);
pde.g_Dn = pde.g_Dn(pt_bot);

pde.exactp = [];
pde.exactux = [];
pde.exactuy = [];
pde.exactu = [];
pde.g_D = [];

% Set initial m0

fprintf('Variable m0\n');
pde.g_R = 1 + mp('0.1') * cos(2 * xbot * mp('pi') + 0.1 * pi);
pde.g_R = linearize_bot(pde.g_R);
m0 = pde.g_R;


pde.g_N = linearize_top(pde.g_N);
pde.g_R = linearize_bot(pde.g_R);
pde.g_RN = (pde.g_RN);
pde.g_Dn = (pde.g_Dn);


[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
uh = soln.u;
ph = soln.p;
u_obs = [uh(IUxTop), uh(IUyTop)];   % Observation corresbone to m0


plt1 = subplot(plot_m, plot_n, 1);
trisurf(elem, node(:,1), node(:,2), uh(IUxNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('sol u', 'FontSize', 14);
view(2);
stokes_plot_solution; % plot the solution


%%  Initial m
m = m0 + mp('0.1') * (sin(xbot * mp('pi') * two) + 0.25);
m = linearize_bot(m);

pde.g_R = m;
plt3 = subplot(plot_m, plot_n, 3);
plot([sft(xbot); sft(xbot)+1], [sft(m0); sft(m0)], '-', ...
     [sft(xbot); sft(xbot)+1], [sft(m); sft(m)], '-x');
legend('m0', 'm');
title('m', 'FontSize', 14);

%% Solve du
% solve L(u, m) = f.
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
um = soln.u;

% plot the difference
plt2 = subplot(plot_m, plot_n, 2);
trisurf(elem, node(:,1), node(:,2), um(IUxNode) - uh(IUxNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('du', 'FontSize', 14);
view(2);

%% Objective
dXidu = two*([um(IUxTop) - u_obs(:,1), xtop*0]);
plt4 = subplot(plot_m, plot_n, 4);
size(xbot);
size(dXidu(:,1));
plot(sft(xbot), sft(dXidu(:,1)), '-');
xlabel('x');
legend('dXidu x');
title('dXidu', 'FontSize', 14);


%% dXidm Adjoint
% Solve adj 
%
pde_adj = pde;
pde_adj.f = 0;
pde_adj.fp = 0;
pde_adj.g_N = -[linearize_top(dXidu(:,1)), xtop*0];
pde_adj.g_R = m;
pde_adj.g_RN = [xbot*0, xbot*0];
pde_adj.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj,option);
ustar = soln.u;

% dXidm = <ustar, dLdm dm>
% dLdm 
dLdm = [um(IUxBot), um(IUyBot)];
%dXidm = h*ustar(IUxBot) .* dLdm;
dXidm = zeros(Nm, 1);
for ii = 1:Nm
    m1 = extend_mid(EI(:, ii));
    dXidm(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                  [ustar(IUxBot), ustar(IUyBot)], ...
                                  dLdm, ...
                                  [m1, m1], option);
end
dXidm = extend_mid(dXidm);

plt5 = subplot(plot_m, plot_n, 5);
plot(sft(xbot), sft(dXidm), '-bx');
title('dXidm', 'FontSize', 14);


%% Iterate to minimize Xi
figure(2);

for k = 1 : max_iteration
    % 
    %  1. Solve u
    %  
    pde_test = pde;
    pde_test.g_R = m;
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
    um = soln.u;
    eqn0 = eqn;

    fprintf('\nIter: %d, Xi: %e, err m: %e, err obs: %e\n', ...
            k, integral_neumann_P2(node, elem, bdFlag, um(IUxTop)-u_obs(:,1), 'repeat', [], option), ...
            norm(m - m0, Inf), norm(um(IUxTop) - u_obs(:,1), Inf) );
    
    % ==================Use Adjoint, with CG==========================
    dXidm = zeros(Nm,1);

    % Primary adj: 
    %    (dLdu)^t u^s = - dXidu^t 
    dXidu = two * [um(IUxTop) - u_obs(:,1), xtop*0];
    pde_adj = pde;
    pde_adj.f = 0;
    pde_adj.fp = 0;
    pde_adj.g_N = -[linearize_top(dXidu(:,1)), xtop*0];
    pde_adj.g_R = linearize_bot(m);
    pde_adj.g_RN = [xbot*0, xbot*0];
    pde_adj.g_Dn = [xbot*0, xbot*0];
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj,option);
    us1 = soln.u;
    
    dLdm = [um(IUxBot) um(IUyBot)];
    for ii = 1:Nm
        m1 = extend_mid(EI(:, ii));
        dXidm(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                      [us1(IUxBot) us1(IUyBot)], ...
                                      dLdm, ...
                                      [m1, m1], ...
                                      option);
    end    
   
    % Stab Matrix
    get_robin_stab_mat;
    
    gamma_stab = 1e-11;

    stokes_info.nx = nx;
    stokes_info.ny = ny;
    stokes_info.xtop = xtop;
    stokes_info.xbot = xbot;
    stokes_info.IUxBot = IUxBot;
    stokes_info.IUyBot = IUyBot;
    stokes_info.IUxTop = IUxTop;
    stokes_info.Mstab = Mstab;
    stokes_info.gamma_stab = gamma_stab;

    dXidm_stab = dXidm;
    
    EI = eye(Nm);
    if true
        % CG iterative
        [dm, flg, relres, niter, resvec] = cgs(@(m1) stokes_hessian(node, elem, stokes_info, bdFlag, m, um, us1, m1, option), ...
                                               dXidm_stab(:), 1e-10, 50);
        fprintf('\tniter: %d, relres: %e\n', niter, relres);
    else
        % inverse
        d2Xidm2 = zeros(Nm);
        for i = 1:Nm
            d2Xidm2(:, i) = stokes_hessian(node, elem, stokes_info, bdFlag, m, um, us1, extend_mid(EI(:,i)), option);
        end
        dm = d2Xidm2 \ dXidm_stab(:);
    end        

    if k <= 4
        subplot(4,1,k)
        plot([xbot(1:n1); xbot(1:n1)+1], [dm; dm], '-bo', ...
             [sft(xbot); sft(xbot)+1], [sft(extend_mid(m-m0)); sft(extend_mid(m-m0))], '-rx');
        legend('delta m', 'm-m0');
        s = sprintf('iter :%d', k);
        title(s);
    end

    % update m
    mbefore = m;
    m = mbefore - extend_mid(dm);
    
end


