# 增量伴随问题的推导

本文解释 [equations.md](equations.md) 中 Gauss--Newton 增量伴随问题的来源，并对应 `LinearStokes/stokes_hessian.m` 的 Hessian-vector product 实现。结构参考 `../FullStokes2d/incremental-adjoint.md`。

## 1. 从伴随问题开始

当前参数 $\beta$ 下，正问题解为 $(\boldsymbol u,p)$，伴随解为 $(\boldsymbol v,r)$。一次伴随问题可抽象写成

$$
\mathcal A(\beta)^*(\boldsymbol v,r)
=
-J_{\boldsymbol u}(\boldsymbol u),
$$

其中 $\mathcal A(\beta)$ 是增量正问题的状态线性化算子。在线性 Stokes 中，体内算子不依赖 $\boldsymbol u$，但边界条件和右端仍依赖当前 $\beta$ 和当前状态。

对本文目标函数

$$
J(\boldsymbol u)
=
\int_{\Gamma_t}
(u_x-u_{{\rm obs},x})^2\,ds,
$$

有

$$
J_{\boldsymbol u}(\boldsymbol u)[\boldsymbol w]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})w_x\,ds.
$$

因此一次伴随顶部牵引为

$$
\frac{\partial\boldsymbol v}{\partial n}-r\boldsymbol n
=
-2(u_x-u_{{\rm obs},x})\boldsymbol e_x.
$$

## 2. 沿参数方向线性化

给定参数方向 $\delta\beta$，考虑

$$
\beta_\epsilon=\beta+\epsilon\delta\beta.
$$

对应的正问题解和伴随解记为

$$
(\boldsymbol u_\epsilon,p_\epsilon),
\qquad
(\boldsymbol v_\epsilon,r_\epsilon).
$$

对 $\epsilon$ 求导得到

$$
(\tilde{\boldsymbol u},\tilde p)
=
\left.
\frac{d}{d\epsilon}
(\boldsymbol u_\epsilon,p_\epsilon)
\right|_{\epsilon=0},
\qquad
(\tilde{\boldsymbol v},\tilde r)
=
\left.
\frac{d}{d\epsilon}
(\boldsymbol v_\epsilon,r_\epsilon)
\right|_{\epsilon=0}.
$$

第一步先解增量正问题。`stokes_hessian.m` 中：

```matlab
pde_du.beta = linearize_bot(m);
pde_du.g_R = [-um(IUxBot), m1, -um(IUyBot), m1];
pde_du.g_D = [xbot*0, xbot*0];
[soln, eqn, info] = StokesP2P1_periodic(node, elem, bdFlag, pde_du, option);
du1 = soln.u;
```

这里 `m1` 是 $\delta\beta$，`du1` 是 $\tilde{\boldsymbol u}$。

## 3. Gauss--Newton 近似

完整 Newton 增量伴随需要线性化一次伴随方程的所有项。线性 Stokes 中，默认 `scheme == 6` 采用 Gauss--Newton 近似：只保留目标函数残差平方项的线性化，不保留与当前伴随变量 $\boldsymbol v$ 相关的底部二阶项。

顶部右端来自

$$
\frac{d}{d\epsilon}
\left[
-2(u_{\epsilon,x}-u_{{\rm obs},x})\boldsymbol e_x
\right]_{\epsilon=0}
=
-2\tilde u_x\boldsymbol e_x.
$$

底部 Robin 右端在 Gauss--Newton 增量伴随中保持零。

## 4. 增量伴随强形式

Gauss--Newton 增量伴随问题为求 $(\tilde{\boldsymbol v},\tilde r)$：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\tilde{\boldsymbol v})+\nabla\tilde r
&=\boldsymbol 0
&&\text{in }\Omega,\\
-\nabla\cdot\tilde{\boldsymbol v}
&=0
&&\text{in }\Omega,\\
\frac{\partial\tilde{\boldsymbol v}}{\partial n}
-\tilde r\boldsymbol n
&=-2\tilde u_x\boldsymbol e_x
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol v}\cdot\boldsymbol n
&=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol v}}{\partial n}
-\tilde r\boldsymbol n
\right)
+\beta\tilde{\boldsymbol v}_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol v},\quad
\left(\frac{\partial\tilde{\boldsymbol v}}{\partial n}
-\tilde r\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

代码中为

```matlab
pde_adj3.f = 0;
pde_adj3.fp = 0;
pde_adj3.g_N = - [2 * linearize_top(du1(IUxTop)), xtop*0];
pde_adj3.beta = linearize_bot(m);
pde_adj3.g_R = [xbot*0, xbot*0];
pde_adj3.g_D = [xbot*0, xbot*0];
[soln, eqn, info] = StokesP2P1_periodic(node, elem, bdFlag, pde_adj3, option);
us3 = soln.u;
```

## 5. Hessian-Vector Product

解出增量伴随速度 `us3` 后，对每个参数基函数 $\phi_i$ 积分：

```matlab
term = integral_robin_P2(node, elem, bdFlag, ...
    [us3(IUxBot), us3(IUyBot)], ...
    [um(IUxBot), um(IUyBot)], ...
    [m2, m2], option);
Hdm(i) = term;
```

对应

$$
(H_{\rm GN}\delta\beta)_i
=
\int_{\Gamma_b}
\phi_i\,
\boldsymbol u_t\cdot\tilde{\boldsymbol v}_t
\,ds.
$$

最后加上很小的稳定化项：

```matlab
Hdm = Hdm + stokes_info.gamma_stab * stokes_info.Mstab * m1(1:n1);
```

默认迭代中 `cgs` 求解

$$
H_{\rm GN}\delta\beta=\nabla J(\beta),
$$

然后

$$
\beta_{\rm new}=\beta-\delta\beta.
$$

## 6. 与完整 Newton 的区别

`scheme == 5` 中保留了完整 Newton 增量伴随的底部源项。一次伴随底部 Robin 条件

$$
\boldsymbol T
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)
+\beta\boldsymbol v_t=0
$$

沿 $\delta\beta$ 线性化时会产生

$$
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol v}^{N}}{\partial n}
-\tilde r^{N}\boldsymbol n
\right)
+\beta\tilde{\boldsymbol v}^{N}_t
=
-\delta\beta\,\boldsymbol v_t.
$$

代码中对应

```matlab
pde_adj2.beta = linearize_bot(beta);
pde_adj2.g_R = - [(v_adj(IUxBot)), beta1, (v_adj(IUyBot)), beta1];
pde_adj2.g_D = [xbot*0, xbot*0];
```

这个分支用于检查完整 Newton Hessian 与有限差分 Hessian，不是默认反演路径。
