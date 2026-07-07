# 线性 Stokes 伴随反演方程与程序对应

本文整理 `../../Stokes/stokes_inversion.m` 中 Robin 系数反演对应的连续方程、
伴随推导和程序实现。这里讨论的是线性 Stokes 方程，不包含非线性黏度或 Glen 流律。

文档中把反演程序里反复出现的四个问题分开写清楚：

1. 正问题：给定当前 Robin 系数 $m$，求状态 $(\boldsymbol u,p)$。
2. 增量正问题：给定参数方向 $\delta m$，求状态方向 $(\tilde{\boldsymbol u},\tilde p)$。
3. 一次伴随问题：由目标函数的一阶状态导数驱动，求 $(\boldsymbol u^*,p^*)$。
4. Gauss--Newton 增量伴随问题：由增量状态的顶部观测导数驱动，求 $(\tilde{\boldsymbol u}^*,\tilde p^*)$。

为突出四个问题的差别，下面只把变化的源项标成红色。

## 1. 基本记号和目标函数

`StokesP2P1_periodic.m` 使用的强形式约定是

$$
-\nabla\cdot(\nabla\boldsymbol u)+\nabla p=\boldsymbol f,
\qquad
-\nabla\cdot\boldsymbol u=f_p
\qquad\text{in }\Omega.
$$

顶部 $\Gamma_t$ 是 Neumann 边界，底部 $\Gamma_b$ 是滑移 Robin 边界，左右边界 $\Gamma_p$ 是周期边界。顶部和底部条件写成

$$
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
=\boldsymbol g_N
\qquad\text{on }\Gamma_t,
$$

$$
\boldsymbol u\cdot\boldsymbol n=g_{Dn},
\qquad
\boldsymbol T
\left(
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
\right)
+
m\boldsymbol u_t
=\boldsymbol g_R
\qquad\text{on }\Gamma_b.
$$

其中 $\boldsymbol T$ 是切向投影，$\boldsymbol u_t=\boldsymbol T\boldsymbol u$。程序内部通过旋转底部自由度实现滑移条件：`pde.g_Dn` 约束法向速度，`pde.g_R` 给出 Robin 系数 $m$，`pde.g_RN` 给出 Robin 边界右端 $\boldsymbol g_R$。

在 `stokes_inversion.m` 中，真实参数记为 `m0`，当前反演参数记为 `m`。程序先用 `m0` 解一次正问题，把顶部速度作为观测：

$$
\boldsymbol u_{\rm obs}
=
\boldsymbol u(m_0)|_{\Gamma_t}.
$$

之后每次迭代给定当前 $m$，解正问题得到 $\boldsymbol u(m)$。目标函数只使用顶部水平速度：

```matlab
dXidu = two*([um(IUxTop) - u_obs(:,1), xtop*0]);
```

对应

$$
\Xi(m)
=
\int_{\Gamma_t}
(u_x-u_{{\rm obs},x})^2\,ds.
$$

这里没有 $\frac12$，因此一阶状态变分是

$$
\delta\Xi(\boldsymbol u)[\boldsymbol w]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})w_x\,ds.
$$

这就是程序中 `dXidu` 前面出现 `two` 的原因。

## 2. 问题一：正问题

给定当前参数 $m$，正问题求 $(\boldsymbol u,p)$：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\boldsymbol u)+\nabla p
&=\color{red}{\boldsymbol f}
&&\text{in }\Omega,\\
-\nabla\cdot\boldsymbol u
&=\color{red}{f_p}
&&\text{in }\Omega,\\
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
&=\color{red}{\boldsymbol g_N}
&&\text{on }\Gamma_t,\\
\boldsymbol u\cdot\boldsymbol n
&=\color{red}{g_{Dn}}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
\right)
+m\boldsymbol u_t
&=\color{red}{\boldsymbol g_R}
&&\text{on }\Gamma_b,\\
\boldsymbol u,\quad
\left(\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

程序中先用 `m0` 生成观测，再用当前 `m` 求反演状态：

```matlab
pde.g_N  = linearize_top(pde.g_N);
pde.g_R  = linearize_bot(pde.g_R);
pde.g_RN = (pde.g_RN);
pde.g_Dn = (pde.g_Dn);

[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
uh = soln.u;
u_obs = [uh(IUxTop), uh(IUyTop)];
```

当前迭代的正问题是

```matlab
pde_test = pde;
pde_test.g_R = m;
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_test,option);
um = soln.u;
```

这里 `um` 对应 $\boldsymbol u(m)$。

## 3. 问题二：增量正问题

设参数方向为 $\delta m$，由它引起的状态方向为 $(\tilde{\boldsymbol u},\tilde p)$。线性 Stokes 的体内算子不依赖状态，所以增量方程的体内源项为零；唯一非零源项来自底部 Robin 项

$$
m\boldsymbol u_t
\quad\Longrightarrow\quad
m\tilde{\boldsymbol u}_t+\delta m\,\boldsymbol u_t.
$$

把 $\delta m\,\boldsymbol u_t$ 移到右端，得到增量正问题：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\tilde{\boldsymbol u})+\nabla\tilde p
&=\color{red}{\boldsymbol 0}
&&\text{in }\Omega,\\
-\nabla\cdot\tilde{\boldsymbol u}
&=\color{red}{0}
&&\text{in }\Omega,\\
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
&=\color{red}{\boldsymbol 0}
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol u}\cdot\boldsymbol n
&=\color{red}{0}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)
+m\tilde{\boldsymbol u}_t
&=\color{red}{-\delta m\,\boldsymbol u_t}
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol u},\quad
\left(\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

这对应 `scheme == 5` 中的增量正问题，也对应 `stokes_hessian.m` 中 Hessian-vector product 的第一步：

```matlab
pde_du.f = 0;
pde_du.fp = 0;
pde_du.g_N = [xtop*0, xtop*0];
pde_du.g_R = linearize_bot(m);
pde_du.g_RN = [-um(IUxBot), m1, -um(IUyBot), m1];
pde_du.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_du,option);
du1 = soln.u;
```

`m1` 是离散方向 $\delta m$。`g_RN` 有四列时，solver 把第 1、2 列相乘作为 $x$ 分量，把第 3、4 列相乘作为 $y$ 分量，因此

$$
[-u_x,m_1,-u_y,m_1]
\quad\Longleftrightarrow\quad
\boldsymbol g_R^{\rm inc}
=
-\delta m\,\boldsymbol u
\quad\text{on }\Gamma_b.
$$

底部法向速度由 `g_Dn=0` 约束，所以这个向量源项在滑移实现中等价于切向右端 $-\delta m\,\boldsymbol u_t$。

## 4. 伴随推导中的配对

为了说明符号来源，把增量正问题的动量和不可压缩方程与伴随变量 $(\boldsymbol v,r)$ 配对：

$$
I
=
\int_\Omega
\left(
-\nabla\cdot(\nabla\tilde{\boldsymbol u})
+\nabla\tilde p
\right)\cdot\boldsymbol v\,dx
-
\int_\Omega
r\,\nabla\cdot\tilde{\boldsymbol u}\,dx.
$$

分部积分并把导数从 $(\tilde{\boldsymbol u},\tilde p)$ 移到 $(\boldsymbol v,r)$ 上，得到

$$
\begin{aligned}
I
&=
\int_\Omega
\left[
-\nabla\cdot(\nabla\boldsymbol v)+\nabla r
\right]\cdot\tilde{\boldsymbol u}\,dx
+
\int_\Omega
\tilde p\,[-\nabla\cdot\boldsymbol v]\,dx\\
&\quad+
\int_{\partial\Omega}
\left[
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)\cdot\tilde{\boldsymbol u}
-
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)\cdot\boldsymbol v
\right]\,ds.
\end{aligned}
$$

因此线性 Stokes 在这个弱形式下的体内伴随算子和原算子相同。伴随边界条件则由上式的边界项和目标函数变分共同决定。

## 5. 问题三：一次伴随问题

一次伴随问题用来消去目标函数一阶变分中的未知状态方向：

$$
\delta\Xi(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})\tilde u_x\,ds.
$$

顶部增量正问题满足齐次 Neumann，因此顶部边界项只剩

$$
\int_{\Gamma_t}
\left(
\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
\right)\cdot\tilde{\boldsymbol u}\,ds.
$$

为了抵消 $\delta\Xi(\boldsymbol u)[\tilde{\boldsymbol u}]$，取一次伴随顶部牵引为

$$
\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
=
\color{red}{-2(u_x-u_{{\rm obs},x})\boldsymbol e_x}
\qquad\text{on }\Gamma_t.
$$

完整的一次伴随问题为求 $(\boldsymbol u^*,p^*)$：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\boldsymbol u^*)+\nabla p^*
&=\color{red}{\boldsymbol 0}
&&\text{in }\Omega,\\
-\nabla\cdot\boldsymbol u^*
&=\color{red}{0}
&&\text{in }\Omega,\\
\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
&=\color{red}{-2(u_x-u_{{\rm obs},x})\boldsymbol e_x}
&&\text{on }\Gamma_t,\\
\boldsymbol u^*\cdot\boldsymbol n
&=\color{red}{0}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
\right)
+m\boldsymbol u^*_t
&=\color{red}{\boldsymbol 0}
&&\text{on }\Gamma_b,\\
\boldsymbol u^*,\quad
\left(\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

程序中对应为：

```matlab
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
```

`us1` 或前面代码里的 `ustar` 就是离散伴随速度 $\boldsymbol u^*$。

## 6. 一阶梯度公式

底部边界项为

$$
-\int_{\Gamma_b}
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)\cdot\boldsymbol u^*\,ds
+
\int_{\Gamma_b}
\left(
\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
\right)\cdot\tilde{\boldsymbol u}\,ds.
$$

在底部不可穿透条件下只需考虑切向分量。把增量 Robin 条件

$$
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)
=
-m\tilde{\boldsymbol u}_t-\delta m\,\boldsymbol u_t
$$

和伴随 Robin 条件

$$
\boldsymbol T
\left(
\frac{\partial\boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
\right)
=
-m\boldsymbol u^*_t
$$

代入，含 $\tilde{\boldsymbol u}_t$ 的项相消，剩下

$$
\delta\Xi(m)[\delta m]
=
\int_{\Gamma_b}
\delta m\,
\boldsymbol u_t\cdot\boldsymbol u^*_t
\,ds.
$$

所以相对于底部 $L^2$ 配对，梯度密度为

$$
g_m
=
\boldsymbol u_t\cdot\boldsymbol u^*_t
\qquad\text{on }\Gamma_b.
$$

程序用每个底部参数基函数 $\phi_i$ 测试这个梯度：

```matlab
dLdm = [um(IUxBot) um(IUyBot)];
for ii = 1:Nm
    m1 = extend_mid(EI(:, ii));
    dXidm(ii) = integral_robin_P2(node, elem, bdFlag, ...
                                  [us1(IUxBot) us1(IUyBot)], ...
                                  dLdm, ...
                                  [m1, m1], ...
                                  option);
end
```

对应

$$
({\rm dXidm})_i
=
\int_{\Gamma_b}
\phi_i\,
\boldsymbol u_t\cdot\boldsymbol u^*_t
\,ds.
$$

`integral_robin_P2` 在 Robin 边界上做 P2 边界积分。第三个函数参数 `[m1,m1]` 表示同一个参数基函数同时乘到速度两个分量上；在滑移实现中，底部法向分量已由 `g_Dn=0` 约束，所以该内积等价于切向内积。

## 7. 有限差分验证

程序随后用中心差分验证伴随梯度：

```matlab
pde_test.g_R = m + ei*deps;
...
pde_test.g_R = m - ei*deps;
...
dXidm_FD(i) = (Xi(m+eps ei)-Xi(m-eps ei))/(2 eps);
```

代码里 `Xi` 通过顶部积分计算：

$$
\Xi(m)
=
\int_{\Gamma_t}
(u_x(m)-u_{{\rm obs},x})^2\,ds.
$$

这部分不参与反演更新，只用于检查伴随梯度符号和大小。

## 8. 问题四：Gauss--Newton 增量伴随问题

在 `scheme == 6` 中，程序用 CGS 求近似 Newton 步：

```matlab
[dm, flg, relres, niter, resvec] = cgs(@(m1) stokes_hessian(..., m1, option), ...
                                       dXidm_stab(:), 1e-10, 50);
```

其中 `stokes_hessian.m` 实现的是 Gauss--Newton Hessian 对方向 $\delta m$ 的作用。给定 $\delta m$，第一步仍是第 3 节的增量正问题，得到 `du1`：

```matlab
pde_du.f = 0;
pde_du.fp = 0;
pde_du.g_N = [xtop*0, xtop*0];
pde_du.g_R = linearize_bot(m);
pde_du.g_RN = [-um(IUxBot), m1, -um(IUyBot), m1];
pde_du.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_du,option);
du1 = soln.u;
```

Gauss--Newton 只保留观测残差的线性化平方项。目标函数二阶状态变分在顶部给出

$$
\delta^2\Xi(\boldsymbol u)
[\tilde{\boldsymbol u},\boldsymbol w]
=
\int_{\Gamma_t}
2\tilde u_x w_x\,ds.
$$

因此增量伴随问题的顶部源项是 $-2\tilde u_x\boldsymbol e_x$，底部 Robin 右端保持零。求 $(\tilde{\boldsymbol u}^*,\tilde p^*)$：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\tilde{\boldsymbol u}^*)+\nabla\tilde p^*
&=\color{red}{\boldsymbol 0}
&&\text{in }\Omega,\\
-\nabla\cdot\tilde{\boldsymbol u}^*
&=\color{red}{0}
&&\text{in }\Omega,\\
\frac{\partial\tilde{\boldsymbol u}^*}{\partial n}
-\tilde p^*\boldsymbol n
&=\color{red}{-2\tilde u_x\boldsymbol e_x}
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol u}^*\cdot\boldsymbol n
&=\color{red}{0}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol u}^*}{\partial n}
-\tilde p^*\boldsymbol n
\right)
+m\tilde{\boldsymbol u}^*_t
&=\color{red}{\boldsymbol 0}
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol u}^*,\quad
\left(\frac{\partial\tilde{\boldsymbol u}^*}{\partial n}
-\tilde p^*\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

程序中对应为：

```matlab
pde_adj3.f = 0;
pde_adj3.fp = 0;
pde_adj3.g_N = - [2 * linearize_top(du1(IUxTop)), xtop*0];
pde_adj3.g_R = linearize_bot(m);
pde_adj3.g_RN = [xbot*0, xbot*0];
pde_adj3.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj3,option);
us3 = soln.u;
```

这里 `us3` 是 Gauss--Newton 增量伴随速度 $\tilde{\boldsymbol u}^*$。然后程序计算

```matlab
term = integral_robin_P2(node, elem, bdFlag, ...
                         [us3(IUxBot), us3(IUyBot)], ...
                         [um(IUxBot), um(IUyBot)], ...
                         [m2, m2], option);
Hdm(i) = term;
```

对应

$$
(H_{\rm GN}\delta m)_i
=
\int_{\Gamma_b}
\phi_i\,
\boldsymbol u_t\cdot\tilde{\boldsymbol u}^*_t
\,ds.
$$

程序最后还加了可选稳定化项：

$$
H_{\rm GN}\delta m
\leftarrow
H_{\rm GN}\delta m
+
\gamma_{\rm stab}M_{\rm stab}\delta m.
$$

对应代码：

```matlab
Hdm = Hdm + stokes_info.gamma_stab * stokes_info.Mstab * m1(1:n1);
```

当前 `stokes_inversion.m` 中设置

```matlab
gamma_stab = 1e-11;
dXidm_stab = dXidm;
```

也就是说右端梯度没有加稳定化梯度，但 Hessian-vector product 中保留了很小的稳定化矩阵项。

## 9. scheme 5 中的完整 Newton 增量伴随

`scheme == 5` 里还有一个用于测试 Hessian 矩阵的完整 Newton 增量伴随问题。它比 Gauss--Newton 增量伴随多了底部源项，这个项来自一次伴随 Robin 条件对参数的线性化：

$$
\boldsymbol T
\left(
\frac{\partial \boldsymbol u^*}{\partial n}
-p^*\boldsymbol n
\right)
+m\boldsymbol u^*_t
=0
\quad\Longrightarrow\quad
\boldsymbol T
\left(
\frac{\partial \tilde{\boldsymbol u}^{*,N}}{\partial n}
-\tilde p^{*,N}\boldsymbol n
\right)
+m\tilde{\boldsymbol u}^{*,N}_t
=
\color{red}{-\delta m\,\boldsymbol u^*_t}.
$$

程序中对应

```matlab
pde_adj2.f = 0;
pde_adj2.fp = 0;
pde_adj2.g_N = - [two * linearize_top(du1(IUxTop)), xtop*0];
pde_adj2.g_R = linearize_bot(m);
pde_adj2.g_RN = - [(us1(IUxBot)), m1, (us1(IUyBot)), m1];
pde_adj2.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde_adj2,option);
us2 = soln.u;
```

因此 `us2` 是完整 Newton 增量伴随，`us3` 是 Gauss--Newton 增量伴随。`scheme == 6` 和 `stokes_hessian.m` 使用的是 `us3`。

完整 Newton 测试中 Hessian 的两个主要项是

$$
\int_{\Gamma_b}
\phi_j\,\boldsymbol u_t\cdot\boldsymbol u^{*,N}_{i,t}\,ds
+
\int_{\Gamma_b}
\phi_j\,\tilde{\boldsymbol u}_{i,t}\cdot\boldsymbol u^*_t\,ds.
$$

代码对应

```matlab
term1 = integral_robin_P2(node, elem, bdFlag, ...
                          [us2(IUxBot), us2(IUyBot)], ...
                          [um(IUxBot), um(IUyBot)], ...
                          [m2, m2], option);

term2 = integral_robin_P2(node, elem, bdFlag, ...
                          [us1(IUxBot), us1(IUyBot)], ...
                          [du1(IUxBot), du1(IUyBot)], ...
                          [m2, m2], option);
```

## 10. 更新公式

CGS 求得 `dm` 后，程序更新

```matlab
m = mbefore - extend_mid(dm);
```

由于 CGS 解的是

$$
H_{\rm GN}\,dm
=
\nabla\Xi(m),
$$

所以更新是

$$
m_{\rm new}
=
m-dm,
$$

即沿 Gauss--Newton 下降方向更新 Robin 系数。

## 11. 程序变量与理论对象对应

| 程序变量 | 理论对象 |
| --- | --- |
| `m0` | 真实 Robin 系数 $m_0$ |
| `m` | 当前反演参数 $m$ |
| `pde.g_R` | Robin 系数 $m$ |
| `pde.g_RN` | Robin 边界右端 $\boldsymbol g_R$ |
| `pde.g_Dn` | 底部法向速度约束，滑移时通常取零 |
| `um` | 当前正问题速度 $\boldsymbol u(m)$ |
| `u_obs` | 顶部观测速度 $\boldsymbol u(m_0)|_{\Gamma_t}$ |
| `dXidu` | 目标函数对状态的导数 $2(u_x-u_{{\rm obs},x})\boldsymbol e_x$ |
| `ustar` / `us1` | 一次伴随速度 $\boldsymbol u^*$ |
| `dXidm` | 梯度 $\boldsymbol u_t\cdot\boldsymbol u^*_t$ 在底部参数基下的系数 |
| `m1` | 参数方向 $\delta m$ |
| `du1` | 增量正问题速度 $\tilde{\boldsymbol u}$ |
| `us2` | 完整 Newton 增量伴随速度 $\tilde{\boldsymbol u}^{*,N}$ |
| `us3` | Gauss--Newton 增量伴随速度 $\tilde{\boldsymbol u}^*$ |
| `stokes_hessian(...)` | Gauss--Newton Hessian-vector product $H_{\rm GN}\delta m$ |
| `integral_neumann_P2` | 顶部 $\Gamma_t$ 上的 P2 边界积分 |
| `integral_robin_P2` | 底部 $\Gamma_b$ 上的 P2 边界积分 |
| `extend_mid` | 把节点参数扩展到 P2 边界节点和边中点 |
| `linearize_top`, `linearize_bot` | 按 solver 需要整理顶部/底部边界 P2 数据 |

## 12. 和非线性 Stokes 文档的关系

本文和 `../nonlinearStokes/theory/adjoint.md` 的结构一致，但线性 Stokes 有两个简化：

1. 体内算子不依赖当前速度，因此没有黏度切线张量 $\mathbb C(\boldsymbol u)$，也没有应力算子的二阶变分。
2. Robin 滑移项是 $m\boldsymbol u_t$，所以参数变分只产生底部项 $\delta m\,\boldsymbol u_t$，梯度自然是底部切向内积 $\boldsymbol u_t\cdot\boldsymbol u^*_t$。

因此，`stokes_inversion.m` 的核心算法可以概括为：

1. 解正问题得到 $\boldsymbol u(m)$；
2. 用顶部速度误差解一次伴随问题得到 $\boldsymbol u^*$；
3. 在底部计算梯度 $\boldsymbol u_t\cdot\boldsymbol u^*_t$；
4. 对每个 CG 方向，解一次增量正问题和一次 Gauss--Newton 增量伴随问题，得到 $H_{\rm GN}\delta m$；
5. 用 $m\leftarrow m-dm$ 更新 Robin 系数。
