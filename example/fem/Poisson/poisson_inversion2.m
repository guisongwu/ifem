%% Poisson inversion

% FIXME: intg on neumann

clear variables;
fprintf('\n\n\n======\nPoisson inversion robin beta\n');
global slope;                           % slab slope

%% Control parameters
slope = 0.1;
h = 0.1;
scheme = 6;
max_iteration = 6;
Xirecord = zeros(max_iteration,1);


% plot dim
plot_m = 2;
plot_n = 3;
option.verb = 0;
option.gNquadorder = 4;
option.use_newton = false;
option.periodic = true;

%% Constants
four = mp('4');
three = mp('3');
two = mp('2');
one = mp('1');

%% Setup domain

[node,elem] = squaremesh([0 1 0 1], h);
bdFlag = setboundary(node,elem, 'Neumann','y==1', 'Robin', 'y==0');
nx = length([0:h:1]);
ny = nx;
n1 = nx-1;

X = reshape(node(:,1), ny, nx);
Y = reshape(node(:,2), ny, nx);
assert(norm(X(1,:)' - [0:h:1]') < 1e-12);
%assert(norm(Y(:,1) - [0:h:1-h]') < 1e-12);
assert(norm(Y(:,1) - [0:h:1]') < 1e-12);

if slope ~= 0  
    % Slab test remap coord
    % Note: this must be done after bdFlag is set     
    fprintf(2, 'Slab test with slope = %f\n', slope);
    node(:,2) = node(:,2) - slope * node(:,1);
end


%% Setup pde and solve 
pde = poisson_data_P2_slab;
warning('\nChange function to data !!!\n');


% bottom
xbot = [0:h:1-h]';              % periodic without last one
ybot = 0 - slope * xbot;
pbot = [xbot ybot];
% top
xtop = xbot;
ytop = 1-h - slope * xbot;
ptop = [xbot ytop];

% top-Neumann
pde.g_N = pde.g_N(ptop);
% bot-Robin
pde.g_RN = pde.g_RN(pbot);

% Set initial m0
if 0
    fprintf('Const m0\n');
    pde.g_R = pde.g_R(pbot); 
else
    fprintf('Variable m0\n');
    pde.g_R = 1 + mp('0.1') * cos(2 * xbot * mp('pi'));
end
m0 = pde.g_R;


[soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde,option);
uh = soln.u;
Uh = reshape(uh, ny, nx); 
u_obs = Uh(ny, 1:nx-1)';   % Observation corresbone to m0


figure;
trisurf(elem, node(:,1), node(:,2), uh(:), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
title('u_0')
colorbar;
view(2);
print('uobsofpoisson','-depsc','-painters');



% if 0
%     int_val = integral_neumann(node, elem, bdFlag, xtop(1:n1).^2.*(1-xtop(1:n1)), xtop(1:n1)*0+1, option);
%     fprintf('test int %e %e\n', int_val, abs(int_val - 1/12));
%     error('stop')
% elseif 0
%     int_val = integral_robin(node, elem, bdFlag, xbot(1:n1).^2.*(1-xbot(1:n1)), xbot(1:n1)*0+1, option);
%     fprintf('test int %e %e\n', int_val, abs(int_val - 1/12));
%     error('stop')
% end

%%  Initial m
%
if 0
    m = m0 + mp('0.02');
elseif 0
    m = m0 + mp('0.1') * sin(xbot * mp('pi'));
elseif 1
    % mo ren this one
    m = m0 + mp('0.1') * (sin(xbot * mp('pi') * two) + 0.25);
elseif 0
    %m = m0 + mp('0.1') * sin(xbot * mp('pi') * three);
    m = m0 + mp('0.05') * sin(xbot * mp('pi') * four);
elseif 0
    m = m0 + mp('0.01') * sin(xbot * mp('pi') * two) ...
        + mp('0.1') * sin(xbot * mp('pi'));
elseif 0
    m = m0 + mp('0.001') * rand(n,1) ...
        + mp('0.0001') * sin(xbot * mp('pi'));
else
    m = m0 + 0.01;
end

pde.g_R = m;

plot([xbot; xbot(end)+h+xbot], [m0; m0], '-o', ...
     [xbot; xbot(end)+h+xbot], [m; m], '--');
legend('m_0', 'm');
title('m_0 and m');
print('mofpoisson','-depsc','-painters');


%% Solve du
% solve L(u, m) = f.
[soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde,option);
um = soln.u;
Um = reshape(um, ny, nx); 

% plot the difference

trisurf(elem, node(:,1), node(:,2), um(:) - uh(:), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
view(2);
print('diffofpoisson','-depsc','-painters');

%% Objective
I = [ny:ny:ny*n1]';      % top without last component (repeated)

dXidu = two*(um(I) - u_obs);
% plt4 = subplot(plot_m, plot_n, 4);
% plot(xbot(1:n1), dXidu, '-o');
% xlabel('x');
% title('dXidu');

%% dXidm Adjoint
% Solve adj 
pde_adj = pde;
pde_adj.f = 0;
pde_adj.g_N = -dXidu;
pde_adj.g_R = m;
pde_adj.g_RN = xbot*0;
[soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj,option);
ustar = soln.u;
Ustar = reshape(ustar, ny, nx);

% dXidm = <ustar, dLdm dm>
% dLdm 
dLdm = Um(1,[1:n1])';
% dXidm = h*Ustar(1,[1:n1])' .* dLdm;
dXidm = zeros(n1,1);
for ii = 1:n1
    m1 = zeros(n1, 1);
    m1(ii) = 1.;
    dXidm(ii) = integral_robin(node, elem, bdFlag, Ustar(1, 1:n1)', dLdm.*m1, option);
end

% Check with 1D
% fprintf('dXidm using adjoint\n');
% fprintf('\t%.8e\n', dXidm);

% plt5 = subplot(plot_m, plot_n, 5);
% plot(xbot, dXidm, '-bx');
% title('dchidm');



%% dXidm FD

deps =  mp(1e-5);

if true
    dXidm_FD = zeros(n1,1);
    for i = 1:n1
        ei = zeros(n1,1);
        ei(i) = 1;

        pde_test = pde;
        pde_test.g_R = m + ei*deps;
        [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
        um1 = soln.u;

        pde_test.g_R = m - ei*deps;
        [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
        um2 = soln.u;
        
        dXidm_FD(i) = sum( (um1(I)-u_obs).^2 - (um2(I)-u_obs).^2 ) *h / deps / 2;
        % fprintf('dXidm(%d) using FD, %.8e\n', i, dXidm_FD(i));
    end
    
%     plt5 = subplot(plot_m, plot_n, 5);
%     hold on;
%     plot(xbot, [dXidm_FD], '-ro');
%     legend('dchidm Adjoint', 'dchidm FD');
    %title('dXidm');
end


%%   Hessian(ii,:) FD
%
Hessian_FD = zeros(n1,n1);
% compute the ii-th row of Hessian
ii = 1;
m1 = zeros(n1,1);
m1(ii) = 1;
for jj = 1 : n1
    m2 = zeros(n1,1);
    m2(jj) = 1;

    pde_test = pde;
    pde_test.g_R = m - m1*deps - m2 *deps;
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
    um00 = soln.u;

    pde_test.g_R = m - m1*deps + m2 *deps;
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
    um01 = soln.u;

    pde_test.g_R = m + m1*deps - m2 *deps;
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
    um10 = soln.u;

    pde_test.g_R = m + m1*deps + m2 *deps;
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
    um11 = soln.u;

    Hessian_FD(ii, jj) = sum( (um11(I)-u_obs).^2 - (um10(I)-u_obs).^2  ...
        - (um01(I)-u_obs).^2 + (um00(I)-u_obs).^2 ) *h / deps/deps/2/2;
end



%% approx Hessian(ii,:) Adjoint
%
Hessian_adj = zeros(n1,n1);
% Incremental:
% (dLdu) du = - dLdm m1
dLdm1 = Um(1, 1:n1)' .* m1;
pde_du = pde;
pde_du.f = 0;
pde_du.g_N = xtop*0;
pde_du.g_R = m;
pde_du.g_RN = -dLdm1;
[soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_du,option);
du1 = soln.u;
Du1 = reshape(du1, ny, nx);

% Second adjoint equation 
pde_adj2 = pde;
pde_adj2.f = 0;
pde_adj2.g_N = - two * Du1(ny, 1:n1)';
pde_adj2.g_R = m;
pde_adj2.g_RN = xbot*0;
[soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj2,option);
us3 = soln.u;
Us3 = reshape(us3, ny, nx);

for jj = 1:n1
    m2 = zeros(n1,1);
    m2(jj) = 1;
    dLdm2 = Um(1, 1:n1)' .* m2; 
    Hessian_adj(ii, jj) = integral_robin(node, elem, bdFlag, Us3(1,1:n1)', dLdm2,option);
end

resultHessian = zeros(3,n1);
resultHessian(1,:) = Hessian_FD(ii,:);
resultHessian(2,:) = Hessian_adj(ii,:);
resultHessian(3,:) = Hessian_FD(ii,:) - Hessian_adj(ii,:);
fprintf('\t no. \t ii-th row of H_diff \t ii-th row of H_adj \t\t misfit\n');
fprintf('\t %3d \t %.8e \t\t %.8e \t %.8e\n', [1:n1; resultHessian]);


% plt6 = subplot(plot_m, plot_n, 6);
% hold on;
% plot(xbot,Hessian_adj(ii,:),'ro-');
% plot(xbot,Hessian_FD(ii,:),'bx--');
% legend('Hessian Adj','Hessian FD');
% hold off;



%% Iterate to minimize Xi
figure(2);

for k = 1 : max_iteration

    %  1. Solve u
    pde_test = pde;
    pde_test.g_R = m;
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
    um = soln.u;


    fprintf('\nIter: %d, Xi: %e, err: %e\n', ...
            k, sum((um(I) - u_obs) .^ 2 * h), norm(m - m0, Inf) );
    Um = reshape(um, ny, nx); 
    Xirecord(k) = sum((um(I) - u_obs) .^ 2 * h);



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
        deps =  mp('1e-9');
        dXidm = zeros(n1,1);
        J = zeros(n1,n1);
        F = zeros(n1,1);

       
        % F
        F = 2 * sqrt(h) * (um(I) - u_obs);
        
        % J = dFdm with FD
        for i = 1:n1
            ei = zeros(n1,1);
            ei(i) = 1;
            
            pde_test = pde;
            pde_test.g_R = m + ei*deps;
            [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
            um1 = soln.u;

            pde_test.g_R = m - ei*deps;
            [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
            um2 = soln.u;
            
            dXidm(i) = sum( (um1(I)-u_obs).^2 - (um2(I)-u_obs).^2 ) *h / deps / 2;
            J(:,i) =  ( um1(I) - um2(I) )' * sqrt(h) / deps / 2;
        end
        
        %disp(det(J'*J));

        if k == 2
            lambda = 1e-9;
        elseif k == 3
            lambda = 1e-10;
        else
            lambda = 1e-10;
        end

        fprintf('\tlambda: %e\n', lambda);

        if false
            dXidm
            2*J'*J
            error('Gauss-Newton stop')
        end

        dm = (J'*J + lambda * eye(n1)) \ dXidm;
        %[dm] = cgs(J'*J + lambda * eye(n), dXidm, 1e-6, 8);
        dm = 1/2*dm;


    elseif scheme == 5

        % --------------------------------------------------------------------------------
        % 
        % Scheme 5: Use Adjoint, only test H_ij
        % 
        % --------------------------------------------------------------------------------
        
        dXidm = zeros(n1,1);
        dXidm_FD = zeros(n1,1);
        d2Xidm2 = zeros(n1,n1);
        d2Xidm2_FD = zeros(n1,n1);

        test_newton = false;
        if test_newton
            warning('Test Newton');
        else
            warning('Test Gauss Newton');
        end
        
        % 
        % Primary adj: 
        %    (dLdu)^t u^s = - dXidu^t 
        % 
        dXidu = two * (um(I) - u_obs);
        pde_adj = pde;
        pde_adj.f = 0;
        pde_adj.g_N = -dXidu;
        pde_adj.g_R = m;
        pde_adj.g_RN = xbot*0;
        [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj,option);
        us1 = soln.u;
        Us1 = reshape(us1, ny, nx);

        dLdm = Um(1, 1:n1)';
        %dXidm = h*Us1(1, 1:n1)' .* dLdm;
        for ii = 1:n1
            m1 = zeros(n1, 1);
            m1(ii) = 1.;
            dXidm(ii) = integral_robin(node, elem, bdFlag, Us1(1, 1:n1)', dLdm.*m1, option);
        end

        for ii = 1:n1
            
            % choose first diriction
            m1 = zeros(n1, 1);
            m1(ii) = 1.;
            
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


            
            if false
                % check du1
                deps =  mp('1e-8');
                
                %u11 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m + m1*deps);
                %u12 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m - m1*deps);

                pde_test = pde;
                pde_test.g_R = m + m1*deps;
                [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
                u11 = soln.u;

                pde_test.g_R = m - m1*deps;
                [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
                u12 = soln.u;
                
                du1_FD = (u11 - u12) / 2 / deps;
                Du1_FD = reshape(du1_FD, ny, nx);

                fprintf('%12.8e %12.8e\n', [Du1(ny, :); ...
                                     Du1_FD(ny, :)]);
                %plot(x, Du1(n, :), '-go', x, Du1_FD(n, :), '-bx');
                error('dbg')
            end

            
            % 
            % Incremental adj: 
            %    L^t u^s = - dXidu^t 
            % 
            % [us2] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), ...
            %                              Us1(1,:) .* m1', m);
            %pde_adj2 = pde;
            pde_adj2.f = 0;
            pde_adj2.g_N = - two * Du1(ny, 1:n1)';
            pde_adj2.g_R = m;
            pde_adj2.g_RN = - Us1(1,1:n1)' .* m1; 

            %warning('us2: dbg');
            % - two * Du1(ny, 1:n1)
            % Us1(1, 1:n1)
            % m1'
            
            [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj2,option);
            us2 = soln.u;
            Us2 = reshape(us2, ny, nx);

            % Gauss Newton           
            % [us3] = poisson_robin2D_mpfr(x, Z(:)*0, - two * Du1(n, :), 0*x, m);
            %pde_adj2 = pde;
            pde_adj2.f = 0;
            pde_adj2.g_N = - two * Du1(ny, 1:n1)';
            pde_adj2.g_R = m;
            pde_adj2.g_RN = xbot*0;
            [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj2,option);
            us3 = soln.u;
            Us3 = reshape(us3, ny, nx);
           
            % choose second direction and test H_ij          
            for jj = 1:n1
                m2 = zeros(n1, 1);
                m2(jj) = 1;

                
                % H_{i,j} = < us2, dLdm(u) dm2  >
                %           + < us1, dLdm(du1) dm2  >
                %           + < us1, d2Ldmdm(u) dm2  >
               
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

                % if ii == 1 && jj == 1
                %     fprintf('%e\n', [norm(Us1)
                %                      norm(Um)
                %                      norm(Du1)
                %                      norm(dLdm2)
                %                      norm(Us2)]);
                %     [term1 term2 term]
                % end

                % if ii <= 2 && jj <= 2
                %     fprintf('term %e %e %e\n', [term1 term2 term1 + term2]);
                %     fprintf('\n');
                % end

                if test_newton
                    % Newton
                    d2Xidm2(ii, jj) = term;
                else
                    % Gauss newton
                    d2Xidm2(ii, jj) = integral_robin(node, elem, bdFlag, Us3(1,1:n1)', dLdm2,option);
                end

                deps =  mp('1e-3');

                if true
                 
                    pde_test = pde;
                    pde_test.g_R = m - m1*deps - m2 *deps;
                    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
                    um00 = soln.u;

                    pde_test.g_R = m - m1*deps + m2 *deps;
                    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
                    um01 = soln.u;

                    pde_test.g_R = m + m1*deps - m2 *deps;
                    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
                    um10 = soln.u;

                    pde_test.g_R = m + m1*deps + m2 *deps;
                    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
                    um11 = soln.u;
                    
                    d2Xidm2_FD(ii, jj) = sum( (um11(I)-u_obs).^2 - (um10(I)-u_obs).^2  ...
                                              - (um01(I)-u_obs).^2 + (um00(I)-u_obs).^2 ) *h / deps/deps/2/2;
                end

                
            end

            deps = 1e-6;

            pde_test = pde;
            pde_test.g_R = m - m1*deps; 
            [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
            um1 = soln.u;

            pde_test.g_R = m + m1*deps; 
            [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_test,option);
            um2 = soln.u;
            
            dXidm_FD(ii) = sum( (um2(I)-u_obs).^2 - (um1(I)-u_obs).^2 ) *h / 2/deps;
        end

        ei = eye(n1); 
        d2Xidm2_func = zeros(n1);
        option.use_newton = test_newton;
        for i = 1:n1
            d2Xidm2_func(:,i) = poisson_hessian(node, elem, bdFlag, m, Um, Us1, ei(:,i), option);
        end

        dm      =  (d2Xidm2 + 1e-9 * eye(n1)) \ dXidm(:);
        dm_func =  d2Xidm2_func \ dXidm(:); 
        dm_FD   = (d2Xidm2_FD + 1e-8 * eye(n1)) \ dXidm_FD(:);


        figure(2);
        if test_newton
            
            plot(xbot, dm, 'o-', ...
                 xbot, dm_func, 'x-', ...
                 xbot, dm_FD, 's-','LineWidth', 2, 'MarkerSize', 10);
            legend('dm adj', 'dm func', 'dm FD');
        end

        if ~test_newton
            plot(xbot, dm, 'o-', xbot, dm_func, 'x-', ...
                 xbot, dm_FD, 's-', xbot, m - m0, '*--', ...
                 'LineWidth', 2, 'MarkerSize', 10);
            legend('dm adj', 'dm func', 'dm FD', 'm-m0');
        elseif 0
            plot(xbot, dm, 'o-', ...
                 xbot, dm_FD, 's-', ...
                 xbot, m - m0, '*--', ...
                 'LineWidth', 2, 'MarkerSize', 10);
            legend('dm adj', 'dm FD', 'm - m0');
        elseif ~test_newton
            plot(xbot, dm, 'o-', ...
                 xbot, m - m0, '*--', xbot, 0*xbot, '-', ...
                 'LineWidth', 2, 'MarkerSize', 10);
            legend('dm adj', 'm - m0', 'zero');
        end

        fprintf('res adj: %e\n', norm(d2Xidm2 * dXidm(:)) / norm(m - m0));
        fprintf('res FD : %e\n', norm(d2Xidm2_FD * dXidm_FD) / norm(m - m0));

        error('stop for testing H_ij')

    elseif scheme == 6
        % --------------------------------------------------------------------------------
        % 
        % Use Adjoint, with CG
        % 
        % --------------------------------------------------------------------------------
        dXidm = zeros(n1,1);

        % 
        % Primary adj: 
        %    (dLdu)^t u^s = - dXidu^t 
        % 
        dXidu = two * (um(I) - u_obs);
        pde_adj = pde;
        pde_adj.f = 0;
        pde_adj.g_N = -dXidu;
        pde_adj.g_R = m;
        pde_adj.g_RN = xbot*0;
        [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde_adj,option);
        us1 = soln.u;
        Us1 = reshape(us1, ny, nx);
        
        dLdm = Um(1, 1:n1)';
        for ii = 1:n1
            m1 = zeros(n1, 1);
            m1(ii) = 1.;
            dXidm(ii) = integral_robin(node, elem, bdFlag, Us1(1, 1:n1)', dLdm.*m1, option);
        end

        if true 
            option.use_newton = false;
            fprintf('Gauss Netwon\n');
        else
            option.use_newton = true;
            fprintf('Netwon\n');
        end
        
        EI = eye(n1);
        if true
            % iterative
            [dm, flg, relres, niter, resvec] = cgs(@(m1) poisson_hessian(node, elem, bdFlag, m, Um, Us1, m1, option), ...
                                                   dXidm(:), 1e-7, 5);
            fprintf('\tniter: %d, relres: %e\n', niter, relres);
         
        else
            % inverse
            d2Xidm2 = zeros(n1);
            for i = 1:n1
                d2Xidm2(:, i) = poisson_hessian(node, elem, bdFlag, m, Um, Us1, EI(:,i), option),
            end
            dm = d2Xidm2 \ dXidm(:);
        end

    end             % end scheme


    if k - 1 <= max_iteration
        subplot(max_iteration,1,k)
        plot([xbot; xbot(end)+h+xbot], [dm; dm], '-bo', ...
             [xbot; xbot(end)+h+xbot], [m - m0; m - m0], '-rx');
        legend('delta m', 'm-m0');
    end

    % update m
    mbefore = m;
    m = mbefore - dm;
    
end
figure(3);
plot(1:max_iteration, Xirecord,'ro-');


