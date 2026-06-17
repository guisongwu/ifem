ii = 1;

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

deps = 1e-4;

fprintf('DXidm   : %e\n', dXidm(ii));
fprintf('DXidm_FD: %e\n', dXidm_FD(ii));

dXidm_2 = 2 * integral_neumann_P2(node, elem, bdFlag, ...
                                  (um(IUxTop)-u_obs(:,1)), du1, [], option);
fprintf('DXidm 2 : %e\n', dXidm_2);     % right






error('dbg dXidm')
