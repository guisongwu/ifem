# 线性 Stokes 伴随反演方程与程序对应

本文整理 `../../../Stokes/stokes_inversion.m` 中 Robin 系数反演使用的连续方程、伴随推导和程序实现对应关系。这里讨论的是线性 Stokes 方程，不包含非线性黏度或 Glen 流律。

## 1. 程序中的正问题

`StokesP2P1_periodic.m` 使用的强形式约定是

$$
-\nabla\cdot(\nabla\boldsymbol u)+\nabla p=\boldsymbol f,
\qquad
-\nabla\cdot\boldsymbol u=f_p
\qquad\text{in }\Omega.
$$

程序中的边界分为：

- 顶部 $\Gamma_t$：Neumann 边界；
- 底部 $\Gamma_b$：Robin 边界；
- 左右边界 $\Gamma_p$：周期边界。

顶部 Neumann 条件为

$$
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
=\boldsymbol g_N
\qquad\text{on }\Gamma_t.
$$

底部在 `option.use_slip=true` 时采用滑移形式。连续层面可写为

$$
\boldsymbol u\cdot\boldsymbol n=0,
\qquad
\boldsymbol T
\left(
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
\right)
+
m\boldsymbol u_t
=\boldsymbol g_R
\qquad\text{on }\Gamma_b,
$$

其中 $\boldsymbol T$ 是切向投影，$\boldsymbol u_t=\boldsymbol T\boldsymbol u$。程序内部通过旋转底部自由度实现这一点：`pde.g_Dn` 固定法向速度，`pde.g_R` 给出切向 Robin 系数。反演参数是底部 Robin 系数 $m$。

在 `stokes_inversion.m` 中，真实参数记为 `m0`，当前反演参数记为 `m`。程序先用 `m0` 求一次正问题，并把顶部速度作为观测：

$$
\boldsymbol u_{\rm obs}
=
\boldsymbol u(m_0)|_{\Gamma_t}.
$$

之后每次迭代给定当前 $m$，解正问题得到 $\boldsymbol u(m)$。

## 2. 目标函数

程序只使用顶部水平速度做观测。代码中

```matlab
dXidu = two*([um(IUxTop) - u_obs(:,1), xtop*0]);
```

对应目标函数

$$
\Xi(m)
=
\int_{\Gamma_t}
(u_x-u_{{\rm obs},x})^2\,ds.
$$

这里没有 $\frac12$，所以对状态变量的一阶变分是

$$
\delta\Xi(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})\tilde u_x\,ds.
$$

这就是程序中 `dXidu` 前面出现 `two` 的原因。

## 3. 增量正问题

设参数方向为 $\delta m$，由它引起的状态方向为 $(\tilde{\boldsymbol u},\tilde p)$。因为线性 Stokes 算子本身不依赖 $\boldsymbol u$，增量正问题只来自底部 Robin 项对 $m$ 的变分。

体内方程为

$$
-\nabla\cdot(\nabla\tilde{\boldsymbol u})+\nabla\tilde p
=\boldsymbol 0,
\qquad
-\nabla\cdot\tilde{\boldsymbol u}=0
\qquad\text{in }\Omega.
$$

顶部为齐次 Neumann：

$$
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-
\tilde p\boldsymbol n
=\boldsymbol 0
\qquad\text{on }\Gamma_t.
$$

底部滑移 Robin 条件由

$$
\boldsymbol u\cdot\boldsymbol n=0,
\qquad
\boldsymbol T
\left(
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
\right)
+
m\boldsymbol u_t
=\boldsymbol g_R
$$

线性化得到

$$
\tilde{\boldsymbol u}\cdot\boldsymbol n=0,
\qquad
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-
\tilde p\boldsymbol n
\right)
+
m\tilde{\boldsymbol u}_t
=
-\delta m\,\boldsymbol u_t
\qquad\text{on }\Gamma_b.
$$

这对应程序中的增量正问题设置：

```matlab
pde_du.g_N  = [xtop*0, xtop*0];
pde_du.g_R  = linearize_bot(m);
pde_du.g_RN = [-um(IUxBot), m1, -um(IUyBot), m1];
pde_du.g_Dn = [xbot*0, xbot*0];
```

`g_RN` 有四列时，solver 把第 1、2 列相乘作为 $x$ 分量，把第 3、4 列相乘作为 $y$ 分量。因此上面代码表示

$$
\boldsymbol g_R^{\rm inc}
=
-\delta m\,\boldsymbol u_t
\qquad\text{on }\Gamma_b.
$$

代码里仍以两个速度分量传入 `g_RN`，但由于底部法向速度由 `g_Dn=0` 固定，这等价于只在切向滑移条件中施加右端。

## 4. 伴随问题

引入伴随变量 $(\boldsymbol u^*,p^*)$。由于线性 Stokes 的主算子在这个弱形式下自伴随，伴随问题使用同一个 Stokes solver，只是边界源项改为目标函数变分给出的顶部牵引。

为了抵消

$$
\delta\Xi(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})\tilde u_x\,ds,
$$

取顶部伴随 Neumann 条件

$$
\frac{\partial\boldsymbol u^*}{\partial n}
-
p^*\boldsymbol n
=
-2(u_x-u_{{\rm obs},x})\boldsymbol e_x
\qquad\text{on }\Gamma_t.
$$

体内伴随方程为

$$
-\nabla\cdot(\nabla\boldsymbol u^*)+\nabla p^*
=\boldsymbol 0,
\qquad
-\nabla\cdot\boldsymbol u^*=0
\qquad\text{in }\Omega.
$$

底部伴随滑移条件保持齐次：

$$
\boldsymbol u^*\cdot\boldsymbol n=0,
\qquad
\boldsymbol T
\left(
\frac{\partial\boldsymbol u^*}{\partial n}
-
p^*\boldsymbol n
\right)
+
m\boldsymbol u^*_t
=\boldsymbol 0
\qquad\text{on }\Gamma_b.
$$

程序中对应为：

```matlab
pde_adj.f    = 0;
pde_adj.fp   = 0;
pde_adj.g_N  = -[linearize_top(dXidu(:,1)), xtop*0];
pde_adj.g_R  = linearize_bot(m);
pde_adj.g_RN = [xbot*0, xbot*0];
pde_adj.g_Dn = [xbot*0, xbot*0];
[soln,eqn,info] = StokesP2P1_periodic(...);
us1 = soln.u;
```

这里 `us1` 就是离散伴随速度 $\boldsymbol u^*$。

## 5. 梯度公式

把增量正问题和伴随问题配对并分部积分。顶部项与目标函数变分相消，体内项由伴随方程消去，底部只剩参数项：

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
                                  [m1, m1], option);
end
```

这对应

$$
({\rm dXidm})_i
=
\int_{\Gamma_b}
\phi_i\,
\boldsymbol u_t\cdot\boldsymbol u^*_t
\,ds.
$$

`integral_robin_P2` 在 Robin 边界上做 P2 边界积分。第三个函数参数 `[m1,m1]` 表示同一个参数基函数同时乘到速度两个分量上；在滑移实现中，底部法向分量已由 `g_Dn=0` 约束，所以该内积等价于切向内积。

## 6. 有限差分验证

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

## 7. Gauss--Newton Hessian-vector product

在 `scheme == 6` 中，程序用 CGS 求解近似 Newton 步：

```matlab
[dm, flg, relres, niter, resvec] = cgs(@(m1) stokes_hessian(..., m1, option), ...
                                       dXidm_stab(:), 1e-10, 50);
```

其中 `stokes_hessian.m` 实现的是 Gauss--Newton Hessian 对方向 $\delta m$ 的作用。

给定方向 $\delta m$，先解增量正问题，得到 $\tilde{\boldsymbol u}$：

```matlab
pde_du.g_N  = [xtop*0, xtop*0];
pde_du.g_R  = linearize_bot(m);
pde_du.g_RN = [-um(IUxBot), m1, -um(IUyBot), m1];
...
du1 = soln.u;
```

目标函数二阶变分为

$$
\delta^2\Xi(\boldsymbol u)
[\tilde{\boldsymbol u},\boldsymbol w]
=
\int_{\Gamma_t}
2\tilde u_x w_x\,ds.
$$

因此增量伴随问题顶部源项为

$$
-2\tilde u_x\boldsymbol e_x.
$$

程序中：

```matlab
pde_adj3.g_N  = - [2 * linearize_top(du1(IUxTop)), xtop*0];
pde_adj3.g_R  = linearize_bot(m);
pde_adj3.g_RN = [xbot*0, xbot*0];
...
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

对应双线性形式

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

## 8. 更新公式

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

即沿 Gauss--Newton 下降方向更新。

## 9. 程序变量与理论对象对应

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
| `us1` / `ustar` | 伴随速度 $\boldsymbol u^*$ |
| `dXidm` | 梯度 $\boldsymbol u_t\cdot\boldsymbol u^*_t$ 在底部参数基下的系数 |
| `m1` | 参数方向 $\delta m$ |
| `du1` | 增量正问题速度 $\tilde{\boldsymbol u}$ |
| `us3` | Gauss--Newton 增量伴随速度 $\tilde{\boldsymbol u}^*$ |
| `stokes_hessian(...)` | Hessian-vector product $H_{\rm GN}\delta m$ |
| `integral_neumann_P2` | 顶部 $\Gamma_t$ 上的 P2 边界积分 |
| `integral_robin_P2` | 底部 $\Gamma_b$ 上的 P2 边界积分 |
| `extend_mid` | 把节点参数扩展到 P2 边界节点和边中点 |
| `linearize_top`, `linearize_bot` | 按 solver 需要整理顶部/底部边界 P2 数据 |

## 10. 和非线性 Stokes 文档的关系

本文和 [equations.md](equations.md) 的结构相同，但有两个简化：

1. 线性 Stokes 的体内算子不依赖当前速度，因此没有黏度切线张量 $\mathbb C(\boldsymbol u)$，也没有应力算子的二阶变分。
2. Robin 滑移项是 $m\boldsymbol u_t$，所以参数变分只产生底部项 $\delta m\,\boldsymbol u_t$，梯度自然是底部切向内积 $\boldsymbol u_t\cdot\boldsymbol u^*_t$。

因此，`stokes_inversion.m` 的核心算法可以概括为：

1. 解正问题得到 $\boldsymbol u(m)$；
2. 用顶部速度误差解伴随问题得到 $\boldsymbol u^*$；
3. 在底部计算梯度 $\boldsymbol u_t\cdot\boldsymbol u^*_t$；
4. 对每个 CG 方向，解一次增量正问题和一次增量伴随问题，得到 $H_{\rm GN}\delta m$；
5. 用 $m\leftarrow m-dm$ 更新 Robin 系数。
