%% Check periodic


figure(2);
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
subplot(1,2,1);
showmesh(node,elem);
findnode(node);

AD = eqn.A;

Ap = eqn.Ap;
bp = eqn.bp;


subplot(1,2,2);
spy(Ap, 10);

%pde.g_D();