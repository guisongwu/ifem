# 伴随问题的推导

本文补充说明 [equations.md](equations.md) 中一次伴随问题的来源。线性 Stokes 情形比 `../FullStokes2d/adjoint.md` 简单：体内状态算子不依赖当前速度，因此伴随算子的体内形式和原算子相同；真正需要处理的是顶部目标函数项和底部 Robin 边界项。

## 1. 从增量正问题开始

固定当前正问题解 $(\boldsymbol u,p)$ 和参数 $\beta$。参数方向 $\delta\beta$ 引起的状态方向 $(\tilde{\boldsymbol u},\tilde p)$ 满足

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\tilde{\boldsymbol u})+\nabla\tilde p
&=\boldsymbol 0,\\
-\nabla\cdot\tilde{\boldsymbol u}
&=0
\end{aligned}
\right.
\qquad\text{in }\Omega,
$$

以及底部线性化 Robin 条件

$$
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)
+\beta\tilde{\boldsymbol u}_t
=
-\delta\beta\,\boldsymbol u_t
\qquad\text{on }\Gamma_b.
$$

目标函数的一阶状态变分是

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})\tilde u_x\,ds.
$$

伴随推导的目标是用 $(\boldsymbol v,r)$ 消去这个顶部的 $\tilde{\boldsymbol u}$，最后只留下底部的 $\delta\beta$ 项。

## 2. 固定弱形式配对

把增量正问题的动量方程和不可压缩方程与伴随变量 $(\boldsymbol v,r)$ 配对：

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
\tilde p[-\nabla\cdot\boldsymbol v]\,dx\\
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

所以体内伴随方程为

$$
-\nabla\cdot(\nabla\boldsymbol v)+\nabla r=\boldsymbol 0,
\qquad
-\nabla\cdot\boldsymbol v=0.
$$

## 3. 顶部边界：观测误差成为伴随牵引

顶部增量正问题满足齐次 Neumann，所以顶部边界项只剩

$$
\int_{\Gamma_t}
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)\cdot\tilde{\boldsymbol u}\,ds.
$$

为了抵消目标函数变分

$$
\int_{\Gamma_t}
2(u_x-u_{{\rm obs},x})\tilde u_x\,ds,
$$

取

$$
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
=
-2(u_x-u_{{\rm obs},x})\boldsymbol e_x
\qquad\text{on }\Gamma_t.
$$

这正对应代码：

```matlab
dJdu = two * [ubeta(IUxTop) - u_obs(:,1), xtop*0];
pde_adj.g_N = -[linearize_top(dJdu(:,1)), xtop*0];
```

## 4. 底部边界：Robin 伴随条件和梯度

底部边界项为

$$
-\int_{\Gamma_b}
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)\cdot\boldsymbol v\,ds
+
\int_{\Gamma_b}
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)\cdot\tilde{\boldsymbol u}\,ds.
$$

由于底部法向分量由 `g_D = 0` 约束，只需看切向分量。增量正问题给出

$$
\boldsymbol T
\left(
\frac{\partial\tilde{\boldsymbol u}}{\partial n}
-\tilde p\boldsymbol n
\right)
=
-\beta\tilde{\boldsymbol u}_t-\delta\beta\,\boldsymbol u_t.
$$

若取伴随底部 Robin 条件

$$
\boldsymbol T
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)
=
-\beta\boldsymbol v_t,
$$

则所有含 $\tilde{\boldsymbol u}_t$ 的项相消，剩下

$$
\delta J(\beta)[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta\,\boldsymbol u_t\cdot\boldsymbol v_t\,ds.
$$

因此梯度密度为

$$
g_\beta=\boldsymbol u_t\cdot\boldsymbol v_t.
$$

程序中对应为

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

## 5. 伴随强形式汇总

一次伴随问题为

$$
\left\{
\begin{aligned}
-\nabla\cdot(\nabla\boldsymbol v)+\nabla r
&=\boldsymbol 0
&&\text{in }\Omega,\\
-\nabla\cdot\boldsymbol v
&=0
&&\text{in }\Omega,\\
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
&=-2(u_x-u_{{\rm obs},x})\boldsymbol e_x
&&\text{on }\Gamma_t,\\
\boldsymbol v\cdot\boldsymbol n
&=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\left(
\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n
\right)
+\beta\boldsymbol v_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol v,\quad
\left(\frac{\partial\boldsymbol v}{\partial n}
-r\boldsymbol n\right)
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

代码中构造为

```matlab
pde_adj = pde;
pde_adj.f = 0;
pde_adj.fp = 0;
pde_adj.g_N = -[linearize_top(dJdu(:,1)), xtop*0];
pde_adj.beta = linearize_bot(beta);
pde_adj.g_R = [xbot*0, xbot*0];
pde_adj.g_D = [xbot*0, xbot*0];
[soln, eqn, info] = StokesP2P1_periodic(node, elem, bdFlag, pde_adj, option);
v_adj = soln.u;
```

## 6. 逻辑小结

线性 Stokes 的一次伴随推导可以概括为：

1. 增量正问题给出参数方向如何进入底部 Robin 条件；
2. 分部积分得到体内伴随算子；
3. 顶部伴随牵引用来抵消目标函数状态变分；
4. 底部伴随 Robin 条件用来消去未知状态方向；
5. 剩下的底部项给出梯度 $\boldsymbol u_t\cdot\boldsymbol v_t$。
