function stokes_inversion_all(mode)
%% Unified Stokes inversion for beta
% Modes:
%   slope     - inclined periodic slab with full diagnostic plots
%   figure    - inclined periodic slab with publication figure export
%   rectangle - horizontal rectangle with constant reference beta

valid_modes = {'slope', 'figure', 'rectangle'};
if nargin < 1
    mode = 'slope';
end
mode = validatestring(mode, valid_modes, mfilename, 'mode');

is_figure_mode = strcmp(mode, 'figure');
is_rectangle_mode = strcmp(mode, 'rectangle');
show_diagnostic_plots = ~is_figure_mode;

% close all; 
% clear variables;
% clf(1);
% clf(2);

fprintf('==============Stokes inversion robin beta===============\n');
global slope h dbg_case dbg_on;                           % slab slope

dbg_on = false;

%% Control parameters
if is_rectangle_mode
    slope = 0;
else
    slope = 0.1;
end
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
if is_figure_mode
    figure;
else
    figure(1);
end
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
Nm = n1; % parameter m with n1 dof
EI = eye(Nm);

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

% Set reference Robin parameter m0.
if is_rectangle_mode
    fprintf('Const m0\n');
    pde.g_R = 0*xbot + 1; %pde.g_R(pt_bot); 
    % pde.g_R = pde.g_R(pt_bot); 
else
    fprintf('Variable m0\n');
    pde.g_R = 1 + mp('0.1') * cos(2 * xbot * mp('pi') + 0.1 * pi);
    pde.g_R = linearize_bot(pde.g_R);
end

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


if ~is_figure_mode
    plt1 = subplot(plot_m, plot_n, 1);
end
trisurf(elem, node(:,1), node(:,2), uh(IUxNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('sol u', 'FontSize', 14);
view(2);

if is_figure_mode
    stokes_plot_solution2;
else
    stokes_plot_solution;
end
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

%%  Initial m
%
if 0
    m = m0 + mp('0.02');
elseif 0
    m = m0 + mp('0.1') * sin(xbot * mp('pi'));
elseif 1
    % mo ren zhe ge
    %m = m0 + mp('1e-5') * (sin(xbot * mp('pi') * two) + 0.25);
    m = m0 + mp('0.1') * (sin(xbot * mp('pi') * two) + 0.25);
    %m = m0 + mp('0.001') * (sin(xbot * mp('pi') * two));
elseif 0
    %m = m0 + mp('0.1') * sin(xbot * mp('pi') * three);
    m = m0 + mp('0.05') * sin(xbot * mp('pi') * four);
elseif 0
    m = m0 + mp('0.01') * sin(xbot * mp('pi') * two) ...
        + mp('0.1') * sin(xbot * mp('pi'));
elseif 0
    m = m0 + mp('0.001') * rand(n,1) ...
        + mp('0.0001') * sin(xbot * mp('pi'));
elseif 0
    m = m0 - 0.1;
elseif 0
    % not working
    m = m + 0.01 * (rand(2*Nm,1) - 0.5);
else
    m = xbot*0 + 1;
end

m = linearize_bot(m);

pde.g_R = m;
if is_figure_mode
    figure;
else
    plt3 = subplot(plot_m, plot_n, 3);
end
%plot(xbot, m0, '-', xbot, m, 'x-');
plot([sft(xbot); sft(xbot)+1], [sft(m0); sft(m0)], '-', ...
     [sft(xbot); sft(xbot)+1], [sft(m); sft(m)], '-x');
legend('m0', 'm');
title('m', 'FontSize', 14);
if is_figure_mode
    print('mofstokes','-depsc','-painters');
end


end

%% Solve du
% solve L(u, m) = f.
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
um = soln.u;


% plot the difference
if is_figure_mode
    figure;
else
    plt2 = subplot(plot_m, plot_n, 2);
end
trisurf(elem, node(:,1), node(:,2), um(IUxNode) - uh(IUxNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('du', 'FontSize', 14);
view(2);
if is_figure_mode
    print('diffofstokes','-depsc','-painters');
end




%% Objective
dXidu = two*([um(IUxTop) - u_obs(:,1), xtop*0]);


if show_diagnostic_plots
    plt4 = subplot(plot_m, plot_n, 4);
    plot(sft(xbot), sft(dXidu(:,1)), '-');
    xlabel('x');
    legend('dXidu x');
    title('dXidu', 'FontSize', 14);
end





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

if show_diagnostic_plots
    fprintf('dXidm using adjoint\n');
    fprintf('\t%.8e\n', dXidm);

    plt5 = subplot(plot_m, plot_n, 5);
    plot(sft(xbot), sft(dXidm), '-bx');
    title('dXidm', 'FontSize', 14);
end



%% dXidm FD
%
%
deps =  mp(1e-5);
dXidm_FD = zeros(Nm,1);

if true
    dXidm_FD = zeros(Nm,1);
    for i = 1:Nm
        ei = extend_mid(EI(:, i));

        pde_test = pde;
        pde_test.g_R = m + ei*deps;
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        um1 = soln.u;

        pde_test.g_R = m - ei*deps;
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        um2 = soln.u;
        
        %dXidm_FD(i) = sum( (um1(IUxTop)-u_obs(:,1)).^2 - (um2(IUxTop)-u_obs(:,1)).^2 ) *h / deps / 2;
        dXidm_FD(i) = (integral_neumann_P2(node, elem, bdFlag, ...
                                           linearize_top(um1(IUxTop)-u_obs(:,1)), 'repeat', [], option) ...
                       - integral_neumann_P2(node, elem, bdFlag, ...
                                             linearize_top(um2(IUxTop)-u_obs(:,1)), 'repeat', [], option)) ...
                                             / deps / 2;
        % fprintf('dXidm(%d) using FD, %.8e\n', i, dXidm_FD(i));
    end
    dXidm_FD = extend_mid(dXidm_FD);
    
    if show_diagnostic_plots
        plt5 = subplot(plot_m, plot_n, 5);
        hold on
        plot(sft(xbot), sft(dXidm_FD), '-ro');
        legend('dXidm adjoint', 'dXidm FD');
    end
end




%% Iterate to minimize Xi
figure(2);

for k = 1 : max_iteration

    %
    % 
    %  1. Solve u
    % 
    % 
    pde_test = pde;
    pde_test.g_R = m;
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
    um = soln.u;
    eqn0 = eqn;

    % warning('um: dbg');
    % % pde_test.g_N'
    % % pde_test.g_R'
    % % pde_test.g_RN'
    % fprintf('\n\n');

    fprintf('\nIter: %d, Xi: %e, err m: %e, err obs: %e\n', ...
            k, integral_neumann_P2(node, elem, bdFlag, um(IUxTop)-u_obs(:,1), 'repeat', [], option), ...
            norm(m - m0, Inf), norm(um(IUxTop) - u_obs(:,1), Inf) );
    


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
        % dXidm = zeros(Nm,1);
        % J = zeros(Nm,Nm);
        % F = zeros(Nm,1);


        % MR = zeros(Nm);
        % for i = 1:Nm
        %     for j = 1:Nm
        %         MR(i, j) = integral_neumann_P2(node, elem, bdFlag, EI(:,i), EI(:,j), [], option);
        %     end
        % end
        % sqrtMR = sqrt(MR);
        

       
        % % F
        % % J = dudm with FD
        % for i = 1:Nm
        %     ei = EI(:, i);
            
        %     pde_test = pde;
        %     pde_test.g_R = m + ei*deps;
        %     [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        %     um1 = soln.u;

        %     pde_test.g_R = m - ei*deps;
        %     [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
        %     um2 = soln.u;
            
        %     dXidm(i) = ( integral_neumann_P2(node, elem, bdFlag, um1(IUxTop)-u_obs(:,1), 'repeat', [], option) ...
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
        %     dXidm
        %     2*J'*J
        %     error('Gauss-Newton stop')
        % end

        % dm = (J'*J + lambda * eye(Nm)) \ dXidm;
        % %[dm] = cgs(J'*J + lambda * eye(n), dXidm, 1e-6, 8);
        % dm = 1/2*dm;


    elseif scheme == 5

        % --------------------------------------------------------------------------------
        % 
        % Scheme 5: Use Adjoint, only test H_ij
        % 
        % --------------------------------------------------------------------------------
        
        dXidm = zeros(Nm,1);
        dXidm_FD = zeros(Nm,1);
        d2Xidm2_NW  = zeros(Nm,Nm);     % Newton
        d2Xidm2_GN = zeros(Nm,Nm);      % Gauss - Newton
        d2Xidm2_FD = zeros(Nm,Nm);


        
        % 
        % Primary adj: 
        %    (dLdu)^t u^s = - dXidu^t 
        % 
        dXidu = two * [um(IUxTop) - u_obs(:,1), xtop*0];
        pde_adj = pde;
        pde_adj.f = 0;
        pde_adj.fp = 0;
        pde_adj.g_N = -[linearize_top(dXidu(:,1)) xtop*0];
        pde_adj.g_R = linearize_bot(m);
        pde_adj.g_RN = [xbot*0, xbot*0];
        pde_adj.g_Dn = [xbot*0, xbot*0];
        [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj,option);
        us1 = soln.u;


        dLdm = [um(IUxBot) um(IUyBot)];
        %dXidm = h*Us1(1, 1:Nm)' .* dLdm;
        for ii = 1:Nm
            m1 = extend_mid(EI(:, ii));
            %dXidm(ii) = integral_robin(node, elem, bdFlag, Us1(1, 1:Nm)', dLdm.*m1, option);
            dXidm(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                          [us1(IUxBot), us1(IUyBot)], ...
                                          dLdm, ...
                                          [m1, m1], ...
                                          option);
        end

        % Test me
        for ii = 1:Nm
        % if true
            % if dbg_case == 1
            %     ii = 1
            % else
            %     ii = 5
            % end
            
            % choose first diriction
            m1 = extend_mid(EI(:, ii));
            
            % 
            % Incremental: 
            %    (dLdu) du = - dLdm m1 
            % 
            dLdm1 = [um(IUxBot), um(IUyBot)].*[m1, m1];
            %pde_du = pde;
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

            %dbgDu1;
            
            % if ii == 1
            %     eqn_m1 = eqn;
            % elseif ii == 2
            %     eqn_m2 = eqn;
            %     error('dbg period');
            % end

            
            % 
            % Incremental adj: 
            %    L^t u^s = - dXidu^t 
            % 
            % [us2] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
            %                              Us1(1,:) .* m1', m);
            %pde_adj2 = pde;
            pde_adj2.f = 0;
            pde_adj2.fp = 0;
            pde_adj2.g_N = - [two * linearize_top(du1(IUxTop)), xtop*0];
            pde_adj2.g_R = linearize_bot(m);
            pde_adj2.g_RN = - [(us1(IUxBot)), m1, (us1(IUyBot)), m1];
            pde_adj2.g_Dn = [xbot*0, xbot*0];
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj2,option);
            us2 = soln.u;

            % 
            % Gauss Newton 
            % 
            % 
            % [us3] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
            %                              0*x, m);
            pde_adj3.f = 0;
            pde_adj3.fp = 0;
            pde_adj3.g_N = - [two * linearize_top(du1(IUxTop)), xtop*0];
            pde_adj3.g_R = linearize_bot(m);
            pde_adj3.g_RN = [xbot*0, xbot*0];
            pde_adj3.g_Dn = [xbot*0, xbot*0];
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj3,option);
            us3 = soln.u;
           

            %error('debug')
            
            % 
            % 
            % choose second direction and test H_ij
            % 
            %
            for jj = 1:Nm
                m2 = extend_mid(EI(:,jj));

                % 
                % H_{i,j} = < us2, dLdm(u) dm2  >
                %           + < us1, dLdm(du1) dm2  >
                %           + < us1, d2Ldmdm(u) dm2  >
                %

                % main term
                dLdm2 = [um(IUxBot), um(IUyBot)] .* [m2, m2];        
                % term1 = dot(Us2(1,1:Nm)', dLdm2);
                term1 = integral_robin_P2(node, elem, bdFlag, ...
                                          [us2(IUxBot), us2(IUyBot)], ...
                                          [um(IUxBot), um(IUyBot)], ...
                                          [m2, m2], option);
                % note: us2 . n | gamma_bot = 0
                

                % second term is much smaller
                dLdm2_du1 = [du1(IUxBot), du1(IUyBot)] .* [m2, m2];        
                % term2 = dot(Us1(1,1:Nm)', dLdm2_du1);
                term2 = integral_robin_P2(node, elem, bdFlag, ...
                                          [us1(IUxBot), us1(IUyBot)], ...
                                          [du1(IUxBot), du1(IUyBot)], ...
                                          [m2, m2], option);

                term3 = 0;
                term = term1 + term2 + term3;

                if ii == 1 && jj == 1
                    % fprintf('%e\n', [norm(Us1)
                    %                  norm(Um)
                    %                  norm(Du1)
                    %                  norm(dLdm2)
                    %                  norm(Us2)]);
                    [term1 term2 term]
                end

                % if ii <= 2 && jj <= 2
                %     fprintf('term %e %e %e\n', [term1 term2 term1 + term2]);
                %     fprintf('\n');
                % end

                % Newton
                d2Xidm2_NW(ii, jj) = term;

                % Gauss newton
                term13 = integral_robin_P2(node, elem, bdFlag, ...
                                           [us3(IUxBot), us3(IUyBot)], ...
                                           [um(IUxBot), um(IUyBot)], ...
                                           [m2, m2], option);
                d2Xidm2_GN(ii, jj) = term13;

                deps =  mp('1e-5');

                if true
                    % compare to FD
                    % um00 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m - m1*deps - m2*deps);
                    % um01 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m - m1*deps + m2*deps);
                    % um10 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m + m1*deps - m2*deps);
                    % um11 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m + m1*deps + m2*deps);

                    pde_test = pde;
                    pde_test.g_R = linearize_bot(m - m1*deps - m2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um00 = soln.u;

                    pde_test.g_R = linearize_bot(m - m1*deps + m2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um01 = soln.u;

                    pde_test.g_R = linearize_bot(m + m1*deps - m2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um10 = soln.u;

                    pde_test.g_R = linearize_bot(m + m1*deps + m2 *deps);
                    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
                    um11 = soln.u;
                    
                    % d2Xidm2_FD(ii, jj) = sum( (um11(I)-u_obs).^2 ...
                    %                           - (um10(I)-u_obs).^2  ...
                    %                           - (um01(I)-u_obs).^2 ...
                    %                           + (um00(I)-u_obs).^2 ) *h / deps/deps/2/2;
                    d2Xidm2_FD(ii, jj) =  (integral_neumann_P2(node, elem, bdFlag, linearize_top(um11(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                                           - integral_neumann_P2(node, elem, bdFlag, linearize_top(um10(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                                           - integral_neumann_P2(node, elem, bdFlag, linearize_top(um01(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                                           + integral_neumann_P2(node, elem, bdFlag, linearize_top(um00(IUxTop)-u_obs(:,1)), 'repeat', [], option))...
                                           / deps/deps/2/2;
                end

                
            end

            % um1 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m - m1*deps);
            % um2 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m + m1*deps);
            
            deps = 1e-6;

            pde_test = pde;
            pde_test.g_R = linearize_bot(m - m1*deps); 
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
            um1 = soln.u;

            pde_test.g_R = linearize_bot(m + m1*deps); 
            [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
            um2 = soln.u;
            
            dXidm_FD(ii) = (integral_neumann_P2(node, elem, bdFlag, linearize_top(um2(IUxTop)-u_obs(:,1)), 'repeat', [], option)...
                            - integral_neumann_P2(node, elem, bdFlag, linearize_top(um1(IUxTop)-u_obs(:,1)), 'repeat', [], option))  / 2/deps;
        end

        % 
        % Stab Matrix
        %
        get_robin_stab_mat;
        
        
        fprintf('res NW: %e\n', norm(d2Xidm2_NW * (m(1:Nm) - m0(1:Nm)) - dXidm(:)) / norm(dXidm));
        fprintf('res GN: %e\n', norm(d2Xidm2_GN * (m(1:Nm) - m0(1:Nm)) - dXidm(:)) / norm(dXidm));
        fprintf('res FD : %e\n', norm(d2Xidm2_FD * (m(1:Nm) - m0(1:Nm)) - dXidm_FD) / norm(dXidm_FD));
        
        %dbgHes;
        dbgDm;
        error('stop for testing H_ij')

    elseif scheme == 6
        % --------------------------------------------------------------------------------
        % 
        % Use Adjoint, with CG
        % 
        % --------------------------------------------------------------------------------
        dXidm = zeros(Nm,1);

        % 
        % Primary adj: 
        %    (dLdu)^t u^s = - dXidu^t 
        % 
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

        %dXidm_stab = dXidm + gamma_stab * Mstab * m(1:n1);
        dXidm_stab = dXidm;

        
        EI = eye(Nm);
        if true
            % iterative
            [dm, flg, relres, niter, resvec] = cgs(@(m1) stokes_hessian(node, elem, stokes_info, bdFlag, m, um, us1, m1, option), ...
                                                   dXidm_stab(:), 1e-10, 50);
            fprintf('\tniter: %d, relres: %e\n', niter, relres);
            %resvec
            % relres
        else
            % inverse
            d2Xidm2 = zeros(Nm);
            for i = 1:Nm
                d2Xidm2(:, i) = stokes_hessian(node, elem, stokes_info, bdFlag, m, um, us1, extend_mid(EI(:,i)), option);
            end
            dm = d2Xidm2 \ dXidm_stab(:);
        end

        if 0
            dXidm(:)
            d2Xidm2
            dm
            error('funcH stop');
        end
        
    end                                 % end scheme


    % if k == 1
    %     subplot(4,1,1)
    %     plot([sft(xbot); sft(xbot)+1], [sft(extend_mid(m0)); sft(extend_mid(m0))], '-rx', ...
    %          [sft(xbot); sft(xbot)+1], [sft(extend_mid(m)); sft(extend_mid(m))], '-bx');
    %     legend('m0', 'm');
    %     title('m0');
    % end
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


end
