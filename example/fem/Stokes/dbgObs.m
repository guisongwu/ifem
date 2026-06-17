

figure(3);

pde.g_R = m0;
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
uh = soln.u;
ph = soln.p;
u_obs = [uh(IUxTop), uh(IUyTop)];   % Observation corresbone to m0


pde.g_R = xbot*0 + 2;
[soln2,eqn2,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
uh2 = soln2.u;
ph2 = soln2.p;
u_obs2 = [uh2(IUxTop), uh2(IUyTop)];   % Observation corresbone to m0

subplot(3,1,1);
plot(sft(xbot), sft(m0), '-o', xbot, pde.g_R, '-o', 'LineWidth', 1);


subplot(3,1,2);
plot(sft(xbot), sft(u_obs), '-o', xbot, sft(u_obs2), '-o', 'LineWidth', 1);


subplot(3,1,3);
patch('faces',elem, 'vertices', node, 'facevertexcdata', uh(IUxNode) - uh2(IUxNode), ...
      'facecolor', 'interp', 'edgecolor', 'none');
axis equal; axis tight; colorbar; view(2);
title('u', 'FontSize', 16);

% subplot(3,1,3);
% %plot(sft(xbot), sft(u_obs), '-o', xbot, u_obs2, '-o', 'LineWidth', 1);
% patch('faces',elem, 'vertices', node, 'facevertexcdata', uh2(IUxNode), ...
%       'facecolor', 'interp', 'edgecolor', 'none');
% axis equal; axis tight; colorbar; view(2);
% title('u', 'FontSize', 16);

