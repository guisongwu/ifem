% check du1
deps =  mp('1e-5');

%u11 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m + m1*deps);
%u12 = poisson_robin2D_mpfr(x, F0(:), gtop, gbot, m - m1*deps);

dbg_on = 0; %true;

pde_test = pde;
pde_test.g_R = linearize_bot(m + m1*deps);
%pde_test.g_R = m1;
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
u11 = soln.u;

% if ii == 1
%     eqn_m1 = eqn;
% elseif ii == 2
%     eqn_m2 = eqn;
% end

% % 
% % (A + dA) * (u + du) ~ F
% % 
% res = eqn1.A_period * (eqn0.sol_period + deps * eqn_d.sol_period) - eqn1.F_period;
% [max(res) min(res)]


pde_test.g_R = linearize_bot(m - m1*deps);
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
u12 = soln.u;

du1_FD = (u11 - u12) / 2 / deps;
%du1_FD = (u11 - um) / deps;

% fprintf('dux %12.8e %12.8e\n', [du1(IUxTop) ...
%                                 du1_FD(IUxTop)]');
% [du1(IUxTop) du1_FD(IUxTop) du1(IUxTop)-du1_FD(IUxTop)]
% [du1(IUxBot) du1_FD(IUxBot) du1(IUxBot)-du1_FD(IUxBot)]
%plot(x, Du1(n, :), '-go', x, Du1_FD(n, :), '-bx');

figure(2);
subplot(4,1,1)
plot(sft(xtop), sft(du1(IUxTop)), '-go', ...
     sft(xtop), sft(du1_FD(IUxTop)), '-bx');
legend('top dux', 'top dux_FD');

subplot(4,1,2)
plot(sft(xtop), sft(du1(IUyTop)), '-go', ...
     sft(xtop), sft(du1_FD(IUyTop)), '-bx');
legend('top duy', 'top duy_FD');

fprintf('max err: %e\n', ...
        max(abs(du1(IUxBot) - du1_FD(IUxBot))));

if 1
    subplot(4,1,3)
    plot([sft(xbot); sft(xbot)+1], [sft(du1(IUxBot)); sft(du1(IUxBot))], '-go', ...
         [sft(xbot); sft(xbot)+1], [sft(du1_FD(IUxBot)); sft(du1_FD(IUxBot))], '-bx');
    legend('bot dux', 'bot dux_FD');

    subplot(4,1,4)
    plot([sft(xbot); sft(xbot)+1], [sft(du1(IUyBot)); sft(du1(IUyBot))], '-go', ...
         [sft(xbot); sft(xbot)+1], [sft(du1_FD(IUyBot)); sft(du1_FD(IUyBot))], '-bx');
    legend('bot duy', 'bot duy_FD');
else
    subplot(4,1,3)
    dd = eqn1.sol_period(IUxBot) - eqn0.sol_period(IUxBot);
    plot([sft(xbot); sft(xbot)+1], [sft(dd); sft(dd)], '-bx');
    legend('bot dux', 'bot dux_FD');

    subplot(4,1,4)
    dd = eqn1.sol_period(IUyBot) - eqn0.sol_period(IUyBot);
    plot([sft(xbot); sft(xbot)+1], [sft(dd); sft(dd)], '-bx');
    legend('bot duy', 'bot duy_FD');
end

if 0
    figure(3);
    subplot(1,2,1)
    trisurf(elem, node(:,1), node(:,2), du1_FD(IUxNode), ...
            'FaceColor', 'interp', 'EdgeColor', 'interp');
    axis equal;
    axis tight;
    colorbar;
    view(2);
    
    subplot(1,2,2)
    trisurf(elem, node(:,1), node(:,2), du1_FD(IUyNode), ...
            'FaceColor', 'interp', 'EdgeColor', 'interp');
    axis equal;
    axis tight;
    colorbar;
    view(2);
end

error('du dbg')

