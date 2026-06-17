
allnode = [unode, ones(Nu,1);
           unode, 2*ones(Nu,1);
           node(:,1), node(:,2) + slope*node(:,1), ones(N,1)];

% allnode(:,2) = allnode(:,2) + slope * allnode(:,1); 

% dA = eqn0.A_period(Nu+1:2*Nu, Nu+1:2*Nu) - eqn1.A_period(Nu+1:2*Nu, Nu+1:2*Nu);
% A0 = eqn0.A_period;

% I = find(sum(dA) ~= 0)';
% unode(I, :)
% J = find(sum(abs(A0(Nu+I,:))) > 1e-12)';
% allnode(J, :)
% eqn1.sol_period(J) - eqn0.sol_period(J)

% I = find(abs(res) > 1e-6);

eps = 1e-12;

I0  = find( abs(allnode(:, 1)) < eps);
II = find( abs(allnode(:, 1) - 1) < eps);

In1 = find( abs(allnode(:, 1) - 1 + h/2) < eps);
In2 = find( abs(allnode(:, 1) - 1 + h) < eps);

In3 = find( abs(allnode(:, 1) - 1 + 3/2*h) < eps);
In4 = find( abs(allnode(:, 1) - 1 + 2*h) < eps);

I1 = find( abs(allnode(:, 1) - h/2) < eps);
I2 = find( abs(allnode(:, 1) - h) < eps);

I3 = find( abs(allnode(:, 1) - 3/2*h) < eps);
I4 = find( abs(allnode(:, 1) - 2*h) < eps);

I5 = find( abs(allnode(:, 1) - 5/2*h) < eps);
I6 = find( abs(allnode(:, 1) - 3*h) < eps);

assert(norm(allnode(In2, 2) - allnode(I0, 2)) < eps);
assert(norm(allnode(I2, 2) - allnode(I0, 2)) < eps );
assert(norm(allnode(I4, 2) - allnode(I0, 2)) < eps);
assert(norm(allnode(I6, 2) - allnode(I0, 2)) < eps );


%% Matrix shift

A1 = eqn_m1.A_period([I0(1:3)], [I0]);
% A2 = eqn_m1.A_period([I2], [I2]);
A2 = eqn_m2.A_period([I2(1:3)], [I2]);

norm(full(A1 - A2))


disp('check 1');
A1 = eqn_m1.A_period([I0], [I0]);
% A2 = eqn_m1.A_period([I2], [I2]);
A2 = eqn_m2.A_period([I2], [I2]);

norm(full(A1 - A2))



disp('check 2');
A1 = eqn_m1.A_period([In1; I0; I1], [In2; In1; I0; I1; I2]);
A2 = eqn_m2.A_period([I1;  I2; I3], [I0;  I1;  I2; I3; I4]);

norm(full(A1 - A2))


A1 = eqn_m1.A_period([In2; In1; I0; I1; I2], [In2; In1; I0; I1; I2]);
A2 = eqn_m2.A_period([I0;  I1;  I2; I3; I4], [I0;  I1;  I2; I3; I4]);

norm(full(A1 - A2))

disp('check 3');
Row1 = [In2; In1; I0; I1; I2;];
Col1 = [In4; In3; In2; In1; I0; I1; I2; I3; I4];

Row2 = [I0; I1; I2; I3; I4;];
Col2 = [In2; In1; I0; I1; I2; I3; I4; I5; I6];

A1 = eqn_m1.A_period(Row1, Col1);
A2 = eqn_m2.A_period(Row2, Col2);

norm(full(A1 - A2))

disp('check length');
length(find(abs(eqn_m1.A_period(Row1, :)) > 1e-12))
length(find(abs(eqn_m1.A_period(Row1, Col1)) > 1e-12))

length(find(abs(eqn_m2.A_period(Row2, :)) > 1e-12))
length(find(abs(eqn_m2.A_period(Row2, Col2)) > 1e-12))


%% Constraint
nny = length(I0);
eqn_m1.A_period(II, I0) + speye(nny)    % error
eqn_m1.A_period(II, II) - speye(nny)


%% Periodic test
Psin = stokes_data_sin_period;
yy = [0:0.01:0.5]';
xx = 0*yy + 0;
nodeL = [xx yy];
nodeR = [xx+1 yy - slope];

disp('check b');
norm(Psin.exactu(nodeL) - Psin.exactu(nodeR))
norm(Psin.exactp(nodeL) - Psin.exactp(nodeR))
norm(Psin.f(nodeL) - Psin.f(nodeR))
norm(Psin.fp(nodeL) - Psin.fp(nodeR))
norm(Psin.g_D(nodeL) - Psin.g_D(nodeR))
norm(Psin.g_N(nodeL) - Psin.g_N(nodeR))
norm(Psin.g_RN(nodeL) - Psin.g_RN(nodeR))
norm(Psin.g_Dn(nodeL) - Psin.g_Dn(nodeR))
