# 线性 Stokes 反问题中的四个方程

本文整理 `LinearStokes/stokes_inversion.m` 中线性 Stokes 底部 Robin 系数反演的连续方程。结构参考 `../FullStokes2d/equations.md`，但这里的体内算子是线性的，不包含 Glen 流律黏度线性化。

反演参数是底部 Robin 系数 $\beta$。源项统一用 `g_*` 命名：

| 程序字段 | 数学对象 |
| :--- | :--- |
| `pde.beta` | Robin 系数 $\beta$，反演参数 |
| `pde.g_R` | Robin 边界右端 $\boldsymbol g_R$ |
| `pde.g_D` | 底部法向速度约束 $g_D$ |
| `pde.g_N` | 顶部 Neumann 数据 $\boldsymbol g_N$ |

## 记号约定

无修饰变量表示当前正问题或伴随问题的主变量：

- $(\boldsymbol u,p)$ 表示给定 $\beta$ 后的正问题速度和压力；
- $(\boldsymbol v,r)$ 表示由当前顶部观测误差驱动的一次伴随速度和伴随压力。

带 $\tilde{\ }$ 的变量表示增量问题中的未知量：

- $(\tilde{\boldsymbol u},\tilde p)$ 表示参数方向 $\delta\beta$ 引起的状态增量；
- $(\tilde{\boldsymbol v},\tilde r)$ 表示 Gauss--Newton Hessian-vector product 中的增量伴随变量。

底部切向投影记为

$$
\boldsymbol T=\boldsymbol I-\boldsymbol n\otimes\boldsymbol n,
\qquad
\boldsymbol w_t=\boldsymbol T\boldsymbol w.
$$

## 几何和离散自由度

脚本先在参考矩形

$$
\widehat\Omega=[0,1]\times[0,0.5]
$$

上生成网格，再通过

$$
y\leftarrow y-sx,\qquad s=\texttt{slope}=0.1
$$

得到倾斜 slab。顶边为 $\Gamma_t$，底边为 $\Gamma_b$，左右边界为周期边界 $\Gamma_p$。

速度采用 $P_2$ 自由度，压力采用 $P_1$ 自由度。底部参数只保留周期独立自由度 `Nbeta`，需要用于边界积分时由 `extend_mid` 扩展到 $P_2$ 边界点。

## 正问题

`StokesP2P1_periodic.m` 使用的强形式约定是

$$
-\nabla\cdot(\nabla\boldsymbol u)+\nabla p=\boldsymbol f,
\qquad
-\nabla\cdot\boldsymbol u=f_p
\qquad\text{in }\Omega.
$$

给定当前参数 $\beta$，正问题求 $(\boldsymbol u,p)$：

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
&=\color{red}{g_D}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n
\right)
+\beta\boldsymbol u_t
&=\color{red}{\boldsymbol g_R}
&&\text{on }\Gamma_b,\\
\boldsymbol u,\quad
\left(\frac{\partial\boldsymbol u}{\partial n}-p\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

对应代码中，当前参数通过 `pde.beta` 传入：

```matlab
pde_test = pde;
pde_test.beta = beta;
[soln, eqn, info] = StokesP2P1_periodic(node, elem, bdFlag, pde_test, option);
ubeta = soln.u;
```

## 目标函数

脚本先用真实参数 `beta0` 解一次正问题，生成顶部速度观测：

```matlab
pde.beta = 1 + 0.1 * cos(2 * xbot * pi + 0.1 * pi);
pde.beta = linearize_bot(pde.beta);
beta0 = pde.beta;

[soln, eqn, info] = StokesP2P1_periodic(node, elem, bdFlag, pde, option);
uh = soln.u;
u_obs = [uh(IUxTop), uh(IUyTop)];
```

目标函数只使用顶部水平速度：

$$
J(\beta)
=
\int_{\Gamma_t}
\left(u_x(\beta)-u_{{\rm obs},x}\right)^2\,ds.
$$

因此状态方向 $\boldsymbol w$ 上的一阶变分为

$$
\delta J(\boldsymbol u)[\boldsymbol w]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})w_x\,ds,
$$

代码中写成

```matlab
dJdu = two * [ubeta(IUxTop) - u_obs(:,1), xtop*0];
```

## 增量正问题

设参数方向为 $\delta\beta$，对应状态方向为 $(\tilde{\boldsymbol u},\tilde p)$。线性 Stokes 体内算子不依赖状态，因此增量方程的体内源项为零。底部 Robin 项

$$
\beta\boldsymbol u_t
$$

线性化后为

$$
\beta\tilde{\boldsymbol u}_t+\delta\beta\,\boldsymbol u_t.
$$

把参数方向项移到右端，得到

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
+\beta\tilde{\boldsymbol u}_t
&=\color{red}{-\delta\beta\,\boldsymbol u_t}
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol u},\quad
\left(\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

`stokes_hessian.m` 中对应为

```matlab
pde_du.beta = linearize_bot(m);
pde_du.g_R = [-um(IUxBot), m1, -um(IUyBot), m1];
pde_du.g_D = [xbot*0, xbot*0];
```

四列 `g_R` 表示逐分量乘积，即

$$
[-u_x,\delta\beta,-u_y,\delta\beta]
\Longleftrightarrow
\boldsymbol g_R^{\rm inc}
=-\delta\beta\,\boldsymbol u.
$$

由于 `g_D = 0` 约束法向速度，这等价于切向右端 $-\delta\beta\,\boldsymbol u_t$。

## 伴随问题

一次伴随问题由目标函数的一阶状态导数驱动，求 $(\boldsymbol v,r)$：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\boldsymbol v)+\nabla r
&=\color{red}{\boldsymbol 0}
&&\text{in }\Omega,\\
-\nabla\cdot\boldsymbol v
&=\color{red}{0}
&&\text{in }\Omega,\\
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
&=\color{red}{-2(u_x-u_{{\rm obs},x})\boldsymbol e_x}
&&\text{on }\Gamma_t,\\
\boldsymbol v\cdot\boldsymbol n
&=\color{red}{0}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)
+\beta\boldsymbol v_t
&=\color{red}{\boldsymbol 0}
&&\text{on }\Gamma_b,\\
\boldsymbol v,\quad
\left(\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

伴随推导和梯度公式见 [adjoint.md](adjoint.md)。

## 梯度公式

伴随方法消去状态方向后得到

$$
\delta J(\beta)[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

因此底部 $L^2$ 配对下的梯度密度为

$$
g_\beta=\boldsymbol u_t\cdot\boldsymbol v_t.
$$

离散实现用每个参数基函数 $\phi_i$ 测试：

```matlab
dLdbeta = [ubeta(IUxBot), ubeta(IUyBot)];
for ii = 1:Nbeta
    beta1 = extend_mid(EI(:, ii));
    dJdbeta(ii) = integral_robin_P2(node, elem, bdFlag, ...
        [v_adj(IUxBot), v_adj(IUyBot)], ...
        dLdbeta, ...
        [beta1, beta1], option);
end
```

## 增量伴随问题

Gauss--Newton Hessian-vector product 中，给定参数方向 $\delta\beta$，先解增量正问题得到 $(\tilde{\boldsymbol u},\tilde p)$，再解增量伴随问题 $(\tilde{\boldsymbol v},\tilde r)$：

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\tilde{\boldsymbol v})+\nabla\tilde r
&=\color{red}{\boldsymbol 0}
&&\text{in }\Omega,\\
-\nabla\cdot\tilde{\boldsymbol v}
&=\color{red}{0}
&&\text{in }\Omega,\\
\frac{\partial\tilde{\boldsymbol v}}{\partial n}
-\tilde r\boldsymbol n
&=\color{red}{-2\tilde u_x\boldsymbol e_x}
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol v}\cdot\boldsymbol n
&=\color{red}{0}
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol v}}{\partial n}
-\tilde r\boldsymbol n
\right)
+\beta\tilde{\boldsymbol v}_t
&=\color{red}{\boldsymbol 0}
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol v},\quad
\left(\frac{\partial\tilde{\boldsymbol v}}{\partial n}
-\tilde r\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

增量伴随问题的来源和 `stokes_hessian.m` 的实现见 [incremental-adjoint.md](incremental-adjoint.md)。

## Gauss--Newton Hessian-Vector Product

`scheme == 6` 中默认使用矩阵自由 Gauss--Newton：

```matlab
[dbeta, flg, relres, niter, resvec] = cgs( ...
    @(beta1) stokes_hessian(node, elem, stokes_info, bdFlag, ...
                            beta, ubeta, v_adj, beta1, option), ...
    dJdbeta_stab(:), 1e-10, 50);
```

也就是求解

$$
H_{\rm GN}\delta\beta=\nabla J(\beta),
$$

然后更新

$$
\beta_{\rm new}=\beta-\delta\beta.
$$

## 符号表

| 程序变量 | 理论对象 |
| :--- | :--- |
| `beta0` | 真实 Robin 系数 $\beta_0$ |
| `beta` | 当前反演参数 $\beta$ |
| `pde.beta` | Robin 系数 $\beta$ |
| `pde.g_R` | Robin 边界右端 $\boldsymbol g_R$ |
| `pde.g_D` | 底部法向速度约束 $g_D$ |
| `ubeta` | 当前正问题速度 $\boldsymbol u(\beta)$ |
| `u_obs` | 顶部观测速度 |
| `dJdu` | 顶部目标函数对速度的导数 |
| `v_adj` | 一次伴随速度 $\boldsymbol v$ |
| `dJdbeta` | 底部参数梯度 |
| `du1` | 增量正问题速度 $\tilde{\boldsymbol u}$ |
| `stokes_hessian(...)` | Gauss--Newton Hessian-vector product |
