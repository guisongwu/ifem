%% Stokes inversion for beta 

% close all; 
% clear variables;
% clf(1);
% clf(2);

fprintf('==============Stokes inversion robin beta===============\n');
global slope h dbg_case dbg_on;                           % slab slope

dbg_on = false;

%% Control parameters
slope = 0.1;
% slope = 0;
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

dbg_case = 0;


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
%assert(norm(Y(:,1) - [0:h:1-h]') < 1e-12);
assert(norm(Y(:,1) - [0:h:.5]') < 1e-12);

N = size(node, 1);
[elem2dof,edge,bdDof] = dofP2(elem);
NE = size(edge,1);
Nu = N + NE; % velocity u using P2
Np = N; % pressure p using P1 
Nbeta = n1; % parameter beta with n1 dof
Nm = Nbeta; % legacy helper scripts still read Nm from the caller workspace
EI = eye(Nbeta);

if slope ~= 0
    %  
    % Slab test remap coord
    % Note: this must be done after bdFlag is set
    % 
    fprintf(2, 'Slab test %f\n', slope);
    node(:,2) = node(:,2) - slope * node(:,1);
    % y = y - slope * x;
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
 

%% 
% Setup pde and solve 
% pde = stokes_data_sin_period;
% pde = stokes_data_P2_period;
pde = stokes_data_grav_period;
warning('\nChange function to data !!!\n');

xbot = [0:h:1-h h/2:h:1]';              % periodic with last one
                                      % first node, then edge center
ybot = 0 - slope * xbot;
% ybot = zeros(size(xbot,1),1);
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

% Set initial beta0
if 0
    fprintf('Const beta0\n');
    pde.g_R = 0*xbot + 1; %pde.g_R(pt_bot); 
    % pde.g_R = pde.g_R(pt_bot); 
else
    fprintf('Variable beta0\n');
    pde.g_R = 1 + mp('0.1') * cos(2 * xbot * mp('pi') + 0.1 * pi);
    pde.g_R = linearize_bot(pde.g_R);
end

pde.g_R = linearize_bot(pde.g_R);
beta0 = pde.g_R;


pde.g_N = linearize_top(pde.g_N);
pde.g_R = linearize_bot(pde.g_R);
pde.g_RN = (pde.g_RN);
pde.g_Dn = (pde.g_Dn);


[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
uh = soln.u;
ph = soln.p;
u_obs = [uh(IUxTop), uh(IUyTop)];   % Observation corresbone to beta0


plt1 = subplot(plot_m, plot_n, 1);
trisurf(elem, node(:,1), node(:,2), uh(IUxNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('sol u', 'FontSize', 14);
view(2);

stokes_plot_solution;
%error('Soultion');


if 0
    int_val = integral_neumann_P2(node, elem, bdFlag, xtop.*(1-xtop), xtop*0+1, [], option);
    fprintf('test int %e %e\n', int_val, abs(int_val - 0.1674979270186815));
    error('stop')
elseif 0
    % int_val = integral_robin_P2(node, elem, bdFlag, xbot.^2.*(1-xbot), xbot*0+1, option);
    % fprintf('test int %e %e\n', int_val, abs(int_val - 1/12));
    int_val = integral_robin_P2(node, elem, bdFlag, [xbot.*(1-xbot), xbot*0+1], [xbot*0+1, xbot.*(1-xbot)], option);
    fprintf('test int %e %e\n', int_val, abs(int_val - 1/6 * 2 * norm([1,slope])));
    error('stop')
end

%%  Initial beta
%
if 0
    beta = beta0 + mp('0.02');
elseif 0
    beta = beta0 + mp('0.1') * sin(xbot * mp('pi'));
elseif 1
    % mo ren zhe ge
    %beta = beta0 + mp('1e-5') * (sin(xbot * mp('pi') * two) + 0.25);
    beta = beta0 + mp('0.1') * (sin(xbot * mp('pi') * two) + 0.25);
    %beta = beta0 + mp('0.001') * (sin(xbot * mp('pi') * two));
elseif 0
    %beta = beta0 + mp('0.1') * sin(xbot * mp('pi') * three);
    beta = beta0 + mp('0.05') * sin(xbot * mp('pi') * four);
elseif 0
    beta = beta0 + mp('0.01') * sin(xbot * mp('pi') * two) ...
        + mp('0.1') * sin(xbot * mp('pi'));
elseif 0
    beta = beta0 + mp('0.001') * rand(n,1) ...
        + mp('0.0001') * sin(xbot * mp('pi'));
elseif 0
    beta = beta0 - 0.1;
elseif 0
    % not working
    beta = beta + 0.01 * (rand(2*Nbeta,1) - 0.5);
else
    beta = xbot*0 + 1;
end

beta = linearize_bot(beta);

pde.g_R = beta;
plt3 = subplot(plot_m, plot_n, 3);
%plot(xbot, beta0, '-', xbot, beta, 'x-');
plot([sft(xbot); sft(xbot)+1], [sft(beta0); sft(beta0)], '-', ...
     [sft(xbot); sft(xbot)+1], [sft(beta); sft(beta)], '-x');
legend('beta0', 'beta');
title('beta', 'FontSize', 14);



%% Solve du
% solve L(u, beta) = f.
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
ubeta = soln.u;


% plot the difference
plt2 = subplot(plot_m, plot_n, 2);
trisurf(elem, node(:,1), node(:,2), ubeta(IUxNode) - uh(IUxNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('du', 'FontSize', 14);
view(2);




%% Objective
dJdu = two*([ubeta(IUxTop) - u_obs(:,1), xtop*0]);


plt4 = subplot(plot_m, plot_n, 4);
size(xbot);
size(dJdu(:,1));
plot(sft(xbot), sft(dJdu(:,1)), '-');
%plot(sft(xbot), sft(dJdu(:,2)), '-');
xlabel('x');
legend('dJdu x');
%legend('dJdu y');
%legend('dJdu x', 'dJdu y');
%ylabel('y');
title('dJdu', 'FontSize', 14);





%% dJdbeta Adjoint
% Solve adj 
%
pde_adj = pde;
pde_adj.f = 0;
pde_adj.fp = 0;
pde_adj.g_N = -[linearize_top(dJdu(:,1)), xtop*0];
pde_adj.g_R = beta;
pde_adj.g_RN = [xbot*0, xbot*0];
pde_adj.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj,option);
v_adj = soln.u;


% dJdbeta = <v_adj, dLdbeta dbeta>
% dLdbeta 
dLdbeta = [ubeta(IUxBot), ubeta(IUyBot)];
%dJdbeta = h*v_adj(IUxBot) .* dLdbeta;
dJdbeta = zeros(Nbeta, 1);
for ii = 1:Nbeta
    beta1 = extend_mid(EI(:, ii));
    dJdbeta(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                  [v_adj(IUxBot), v_adj(IUyBot)], ...
                                  dLdbeta, ...
                                  [beta1, beta1], option);
end
dJdbeta = extend_mid(dJdbeta);

% Check with 1D
fprintf('dJdbeta using adjoint\n');
fprintf('\t%.8e\n', dJdbeta);

plt5 = subplot(plot_m, plot_n, 5);
plot(sft(xbot), sft(dJdbeta), '-bx');
title('dJdbeta', 'FontSize', 14);



%% dJdbeta FD
%
%
deps =  mp(1e-5);
dJdbeta_FD = zeros(Nbeta,1);

if true
    dJdbeta_FD = zeros(Nbeta,1);
    for i = 1:Nbeta
        ei = extend_mid(EI(:, i));

        pde_test = pde;
        pde_test.g_R = beta + ei*deps;
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        um1 = soln.u;

        pde_test.g_R = beta - ei*deps;
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        um2 = soln.u;
        
        %dJdbeta_FD(i) = sum( (um1(IUxTop)-u_obs(:,1)).^2 - (um2(IUxTop)-u_obs(:,1)).^2 ) *h / deps / 2;
        dJdbeta_FD(i) = (integral_neumann_P2(node, elem, bdFlag, ...
                                           linearize_top(um1(IUxTop)-u_obs(:,1)), 'repeat', [], option) ...
                       - integral_neumann_P2(node, elem, bdFlag, ...
                                             linearize_top(um2(IUxTop)-u_obs(:,1)), 'repeat', [], option)) ...
                                             / deps / 2;
        % fprintf('dJdbeta(%d) using FD, %.8e\n', i, dJdbeta_FD(i));
    end
    dJdbeta_FD = extend_mid(dJdbeta_FD);
    
    plt5 = subplot(plot_m, plot_n, 5);
    hold on
    plot(sft(xbot), sft(dJdbeta_FD), '-ro');
    %plot(xbot, dJdbeta_FD, '-ro');
    legend('dJdbeta adjoint', 'dJdbeta FD');
    %title('dJdbeta');
end




%% Iterate to minimize J
figure(2);

for k = 1 : max_iteration

    %
    % 
    %  1. Solve u
    % 
    % 
    pde_test = pde;
    pde_test.g_R = beta;
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
    ubeta = soln.u;
    eqn0 = eqn;

    % warning('ubeta: dbg');
    % % pde_test.g_N'
    % % pde_test.g_R'
    % % pde_test.g_RN'
    % fprintf('\n\n');

    fprintf('\nIter: %d, J: %e, err beta: %e, err obs: %e\n', ...
            k, integral_neumann_P2(node, elem, bdFlag, ubeta(IUxTop)-u_obs(:,1), 'repeat', [], option), ...
            norm(beta - beta0, Inf), norm(ubeta(IUxTop) - u_obs(:,1), Inf) );
    


    % 
    % 
    % Scheme: 
    %    4: Gauss newton with Jacobian
    %    5: Newton with FD 
    %    6: Adjoint
    % 
    %
    if scheme == 4
        % --------------------------------------------------------------------------------
        % 
        % Scheme 4: Gauss newton with Jacobian
        % 
        % --------------------------------------------------------------------------------

        % 
        % 
        % Note: not implemented 
        % 
        % 
        
        % deps =  mp('1e-9');
        % dJdbeta = zeros(Nbeta,1);
        % J = zeros(Nbeta,Nbeta);
        % F = zeros(Nbeta,1);


        % MR = zeros(Nbeta);
        % for i = 1:Nbeta
        %     for j = 1:Nbeta
        %         MR(i, j) = integral_neumann_P2(node, elem, bdFlag, EI(:,i), EI(:,j), [], option);
        %     end
        % end
        % sqrtMR = sqrt(MR);
        

       
        % % F
        % % J = dudbeta with FD
        % for i = 1:Nbeta
        %     ei = EI(:, i);
            
        %     pde_test = pde;
        %     pde_test.g_R = beta + ei*deps;
        %     [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        %     um1 = soln.u;

        %     pde_test.g_R = beta - ei*deps;
        %     [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        %     um2 = soln.u;
            
        %     dJdbeta(i) = ( integral_neumann_P2(node, elem, bdFlag, um1(IUxTop)-u_obs(:,1), 'repeat', [], option) ...
        %                  - integral_neumann_P2(node, elem, bdFlag, um2(IUxTop)-u_obs(:,1), 'repeat', [], option)) ...
        %                  / deps / 2;
        %     J(:,i) =  sqrtMR * ( um1(IUxTop) - um2(IUxTop) ) / deps / 2;
        % end
        
        % %disp(det(J'*J));

        % if k == 2
        %     lambda = 1e-9;
        % elseif k == 3
        %     lambda = 1e-10;
        % else
        %     lambda = 1e-10;
        % end

        % fprintf('\tlambda: %e\n', lambda);

        % if false
        %     dJdbeta
        %     2*J'*J
        %     error('Gauss-Newton stop')
        % end

        % dbeta = (J'*J + lambda * eye(Nbeta)) \ dJdbeta;
        % %[dbeta] = cgs(J'*J + lambda * eye(n), dJdbeta, 1e-6, 8);
        % dbeta = 1/2*dbeta;


    elseif scheme == 5

        % --------------------------------------------------------------------------------
        % 
        % Scheme 5: Use Adjoint, only test H_ij
        % 
        % --------------------------------------------------------------------------------
        
        dJdbeta = zeros(Nbeta,1);
        dJdbeta_FD = zeros(Nbeta,1);
        d2Jdbeta2_NW  = zeros(Nbeta,Nbeta);     % Newton
        d2Jdbeta2_GN = zeros(Nbeta,Nbeta);      % Gauss - Newton
        d2Jdbeta2_FD = zeros(Nbeta,Nbeta);


        
        % 
        % Primary adj: 
        %    (dLdu)^t u^s = - dJdu^t 
        % 
        dJdu = two * [ubeta(IUxTop) - u_obs(:,1), xtop*0];
        pde_adj = pde;
        pde_adj.f = 0;
        pde_adj.fp = 0;
        pde_adj.g_N = -[linearize_top(dJdu(:,1)) xtop*0];
        pde_adj.g_R = linearize_bot(beta);
        pde_adj.g_RN = [xbot*0, xbot*0];
        pde_adj.g_Dn = [xbot*0, xbot*0];
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj,option);
        v_adj = soln.u;


        dLdbeta = [ubeta(IUxBot) ubeta(IUyBot)];
        %dJdbeta = h*Us1(1, 1:Nbeta)' .* dLdbeta;
        for ii = 1:Nbeta
            beta1 = extend_mid(EI(:, ii));
            %dJdbeta(ii) = integral_robin(node, elem, bdFlag, Us1(1, 1:Nbeta)', dLdbeta.*beta1, option);
            dJdbeta(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                          [v_adj(IUxBot), v_adj(IUyBot)], ...
                                          dLdbeta, ...
                                          [beta1, beta1], ...
                                          option);
        end

        % Test me
        for ii = 1:Nbeta
        % if true
            % if dbg_case == 1
            %     ii = 1
            % else
            %     ii = 5
            % end
            
            % choose first diriction
            beta1 = extend_mid(EI(:, ii));
            
            % 
            % Incremental: 
            %    (dLdu) du = - dLdbeta beta1 
            % 
            dLdbeta1 = [ubeta(IUxBot), ubeta(IUyBot)].*[beta1, beta1];
            %pde_du = pde;
            pde_du.f = 0;
            pde_du.fp = 0;
            pde_du.g_N = [xtop*0, xtop*0];
            pde_du.g_R = linearize_bot(beta);
            %pde_du.g_RN = [-linearize_bot(dLdbeta1(:,1)), -linearize_bot(dLdbeta1(:,2))];
            pde_du.g_RN = [-ubeta(IUxBot), beta1, -ubeta(IUyBot), beta1]; 
            pde_du.g_Dn = [xbot*0, xbot*0];
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_du,option);
            % [du1] = poisson_robin2D_mpfr(x, Z(:)*0, 0*x, -dLdbeta1, beta);
            du1 = soln.u;

            %dbgDu1;
            
            % if ii == 1
            %     eqn_m1 = eqn;
            % elseif ii == 2
            %     eqn_m2 = eqn;
            %     error('dbg period');
            % end

            
            % 
            % Incremental adj: 
            %    L^t u^s = - dJdu^t 
            % 
            % [v_newton] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
            %                              Us1(1,:) .* beta1', beta);
            %pde_adj2 = pde;
            pde_adj2.f = 0;
            pde_adj2.fp = 0;
            pde_adj2.g_N = - [two * linearize_top(du1(IUxTop)), xtop*0];
            pde_adj2.g_R = linearize_bot(beta);
            pde_adj2.g_RN = - [(v_adj(IUxBot)), beta1, (v_adj(IUyBot)), beta1];
            pde_adj2.g_Dn = [xbot*0, xbot*0];
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj2,option);
            v_newton = soln.u;

            % 
            % Gauss Newton 
            % 
            % 
            % [v_gn] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
            %                              0*x, beta);
            pde_adj3.f = 0;
            pde_adj3.fp = 0;
            pde_adj3.g_N = - [two * linearize_top(du1(IUxTop)), xtop*0];
            pde_adj3.g_R = linearize_bot(beta);
            pde_adj3.g_RN = [xbot*0, xbot*0];
            pde_adj3.g_Dn = [xbot*0, xbot*0];
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj3,option);
            v_gn = soln.u;
           

            %error('debug')
            
            % 
            % 
            % choose second direction and test H_ij
            % 
            %
            for jj = 1:Nbeta
                beta2 = extend_mid(EI(:,jj));

                % 
                % H_{i,j} = < v_newton, dLdbeta(u) dbeta2  >
                %           + < v_adj, dLdbeta(du1) dbeta2  >
                %           + < v_adj, d2Ldbetadbeta(u) dbeta2  >
                %

                % main term
                dLdbeta2 = [ubeta(IUxBot), ubeta(IUyBot)] .* [beta2, beta2];        
                % term1 = dot(Us2(1,1:Nbeta)', dLdbeta2);
                term1 = integral_robin_P2(node, elem, bdFlag, ...
                                          [v_newton(IUxBot), v_newton(IUyBot)], ...
                                          [ubeta(IUxBot), ubeta(IUyBot)], ...
                                          [beta2, beta2], option);
                % note: v_newton . n | gamma_bot = 0
                

                % second term is much smaller
                dLdbeta2_du1 = [du1(IUxBot), du1(IUyBot)] .* [beta2, beta2];        
                % term2 = dot(Us1(1,1:Nbeta)', dLdbeta2_du1);
                term2 = integral_robin_P2(node, elem, bdFlag, ...
                                          [v_adj(IUxBot), v_adj(IUyBot)], ...
                                          [du1(IUxBot), du1(IUyBot)], ...
                                          [beta2, beta2], option);

                term3 = 0;
                term = term1 + term2 + term3;

                if ii == 1 && jj == 1
                    % fprintf('%e\n', [norm(Us1)
                    %                  norm(Um)
                    %                  norm(Du1)
                    %                  norm(dLdbeta2)
                    %                  norm(Us2)]);
                    [term1 term2 term]
                end

                % if ii <= 2 && jj <= 2
                %     fprintf('term %e %e %e\n', [term1 term2 term1 + term2]);
                %     fprintf('\n');
                % end

                % Newton
                d2Jdbeta2_NW(ii, jj) = term;

                % Gauss newton
                term13 = integral_robin_P2(node, elem, bdFlag, ...
                                           [v_gn(IUxBot), v_gn(IUyBot)], ...
                                           [ubeta(IUxBot), ubeta(IUyBot)], ...
                                           [beta2, beta2], option);
                d2Jdbeta2_GN(ii, jj) = term13;

                deps =  mp('1e-5');

                if true
                    % compare to FD
                    % um00 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, beta - beta1*deps - beta2*deps);
                    % um01 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, beta - beta1*deps + beta2*deps);
                    % um10 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, beta + beta1*deps - beta2*deps);
                    % um11 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, beta + beta1*deps + beta2*deps);

                    pde_test = pde;
                    pde_test.g_R = linearize_bot(beta - beta1*deps - beta2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um00 = soln.u;

                    pde_test.g_R = linearize_bot(beta - beta1*deps + beta2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um01 = soln.u;

                    pde_test.g_R = linearize_bot(beta + beta1*deps - beta2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um10 = soln.u;

                    pde_test.g_R = linearize_bot(beta + beta1*deps + beta2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um11 = soln.u;
                    
                    % d2Jdbeta2_FD(ii, jj) = sum( (um11(I)-u_obs).^2 ...
                    %                           - (um10(I)-u_obs).^2  ...
                    %                           - (um01(I)-u_obs).^2 ...
                    %                           + (um00(I)-u_obs).^2 ) *h / deps/deps/2/2;
                    d2Jdbeta2_FD(ii, jj) =  (integral_neumann_P2(node, elem, bdFlag, linearize_top(um11(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                                           - integral_neumann_P2(node, elem, bdFlag, linearize_top(um10(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                                           - integral_neumann_P2(node, elem, bdFlag, linearize_top(um01(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                                           + integral_neumann_P2(node, elem, bdFlag, linearize_top(um00(IUxTop)-u_obs(:,1)), 'repeat', [], option))...
                                           / deps/deps/2/2;
                end

                
            end

            % um1 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, beta - beta1*deps);
            % um2 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, beta + beta1*deps);
            
            deps = 1e-6;

            pde_test = pde;
            pde_test.g_R = linearize_bot(beta - beta1*deps); 
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
            um1 = soln.u;

            pde_test.g_R = linearize_bot(beta + beta1*deps); 
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
            um2 = soln.u;
            
            dJdbeta_FD(ii) = (integral_neumann_P2(node, elem, bdFlag, linearize_top(um2(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                            - integral_neumann_P2(node, elem, bdFlag, linearize_top(um1(IUxTop)-u_obs(:,1)), 'repeat', [], option))  / 2/deps;
        end

        % 
        % Stab Matrix
        %
        get_robin_stab_mat;
        
        
        fprintf('res NW: %e\n', norm(d2Jdbeta2_NW * (beta(1:Nbeta) - beta0(1:Nbeta)) - dJdbeta(:)) / norm(dJdbeta));
        fprintf('res GN: %e\n', norm(d2Jdbeta2_GN * (beta(1:Nbeta) - beta0(1:Nbeta)) - dJdbeta(:)) / norm(dJdbeta));
        fprintf('res FD : %e\n', norm(d2Jdbeta2_FD * (beta(1:Nbeta) - beta0(1:Nbeta)) - dJdbeta_FD) / norm(dJdbeta_FD));
        
        %dbgHes;
        dbgDm;
        error('stop for testing H_ij')

    elseif scheme == 6
        % --------------------------------------------------------------------------------
        % 
        % Use Adjoint, with CG
        % 
        % --------------------------------------------------------------------------------
        dJdbeta = zeros(Nbeta,1);

        % 
        % Primary adj: 
        %    (dLdu)^t u^s = - dJdu^t 
        % 
        dJdu = two * [ubeta(IUxTop) - u_obs(:,1), xtop*0];
        pde_adj = pde;
        pde_adj.f = 0;
        pde_adj.fp = 0;
        pde_adj.g_N = -[linearize_top(dJdu(:,1)), xtop*0];
        pde_adj.g_R = linearize_bot(beta);
        pde_adj.g_RN = [xbot*0, xbot*0];
        pde_adj.g_Dn = [xbot*0, xbot*0];
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj,option);
        v_adj = soln.u;
        
        dLdbeta = [ubeta(IUxBot) ubeta(IUyBot)];
        for ii = 1:Nbeta
            beta1 = extend_mid(EI(:, ii));
            dJdbeta(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                          [v_adj(IUxBot) v_adj(IUyBot)], ...
                                          dLdbeta, ...
                                          [beta1, beta1], ...
                                          option);
        end

        
        % 
        % Stab Matrix
        %
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

        %dJdbeta_stab = dJdbeta + gamma_stab * Mstab * beta(1:n1);
        dJdbeta_stab = dJdbeta;

        
        EI = eye(Nbeta);
        if true
            % iterative
            [dbeta, flg, relres, niter, resvec] = cgs(@(beta1) stokes_hessian(node, elem, stokes_info, bdFlag, beta, ubeta, v_adj, beta1, option), ...
                                                   dJdbeta_stab(:), 1e-10, 50);
            fprintf('\tniter: %d, relres: %e\n', niter, relres);
            %resvec
            % relres
        else
            % inverse
            d2Jdbeta2 = zeros(Nbeta);
            for i = 1:Nbeta
                d2Jdbeta2(:, i) = stokes_hessian(node, elem, stokes_info, bdFlag, beta, ubeta, v_adj, extend_mid(EI(:,i)), option);
            end
            dbeta = d2Jdbeta2 \ dJdbeta_stab(:);
        end

        if 0
            dJdbeta(:)
            d2Jdbeta2
            dbeta
            error('funcH stop');
        end
        
    end                                 % end scheme


    % if k == 1
    %     subplot(4,1,1)
    %     plot([sft(xbot); sft(xbot)+1], [sft(extend_mid(beta0)); sft(extend_mid(beta0))], '-rx', ...
    %          [sft(xbot); sft(xbot)+1], [sft(extend_mid(beta)); sft(extend_mid(beta))], '-bx');
    %     legend('beta0', 'beta');
    %     title('beta0');
    % end
    if k <= 4
        subplot(4,1,k)
        plot([xbot(1:n1); xbot(1:n1)+1], [dbeta; dbeta], '-bo', ...
             [sft(xbot); sft(xbot)+1], [sft(extend_mid(beta-beta0)); sft(extend_mid(beta-beta0))], '-rx');
        legend('delta beta', 'beta-beta0');
        s = sprintf('iter :%d', k);
        title(s);
    end

    % update beta
    betabefore = beta;
    beta = betabefore - extend_mid(dbeta);
    
end


