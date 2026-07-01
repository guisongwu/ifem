# First-order 冰流模型伴随反问题方程

本文整理 `FirstOrderP2.m` 和 `FOAdjInvSlabBed.m` 中使用的
two-dimensional first-order/BP slab 模型及其伴随反演方程。目标是说明：

1. 正问题；
2. 增量正问题；
3. 伴随问题；
4. Gauss--Newton Hessian-vector product。

这里讨论的是连续弱形式和它与 MATLAB 实现的对应关系。`FirstOrderP2.m`
使用 P2 三角元离散水平速度，`FOAdjInvSlabBed.m` 用顶部水平速度观测反演底部
摩擦系数。

## 几何与记号

区域为倾斜周期 slab：

$$
\Omega
=
\{(x,z):0<x<L,\ b(x)<z<s(x)\},
\qquad
b(x)=-\alpha x,\quad s(x)=H-\alpha x.
$$

边界分为

$$
\partial\Omega=\Gamma_t\cup\Gamma_b\cup\Gamma_p,
$$

其中 $\Gamma_t$ 是上表面，$\Gamma_b$ 是底部，$\Gamma_p$ 是左右周期边界。
周期条件使用坐标

$$
\xi=(x,\ z+\alpha x)
$$

识别左右边界，因此右边界点 $(L,z)$ 与左边界点
$(0,z+\alpha L)$ 周期等价。

first-order/BP 横截面模型只求水平速度

$$
u=u(x,z).
$$

诊断竖直速度 $w$ 和压力 $p$ 在 `FirstOrderP2.m` 中后处理得到，不参与反问题的主
方程。底部单位切向量记为

$$
\boldsymbol t=(t_x,t_z)^T,
\qquad
\chi=t_x^2.
$$

由于未知量只有水平速度，底部滑移项在标量弱形式中带有几何因子 $\chi$。

## 黏度和正问题

FO 模型使用 Glen 型有效黏度

$$
\eta(u)
=
\frac12 A^{-1/n}
\left(
\varepsilon_{\rm II}(u)+\varepsilon_{\rm reg}^2
\right)^{\frac{1-n}{2n}},
$$

其中

$$
\varepsilon_{\rm II}(u)
=
u_x^2+\frac14 u_z^2.
$$

记

$$
\tau_x(u)=4\eta(u)u_x,
\qquad
\tau_z(u)=\eta(u)u_z.
$$

正问题强形式为：给定底部摩擦系数 $\beta$，求 $u$，满足

$$
\left\{
\begin{aligned}
-\partial_x\tau_x(u)-\partial_z\tau_z(u) &= f
&&\text{in }\Omega,\\
\tau_x(u)n_x+\tau_z(u)n_z &=0
&&\text{on }\Gamma_t,\\
\tau_x(u)n_x+\tau_z(u)n_z
+\beta |u t_x|^{m-1}u\chi &=0
&&\text{on }\Gamma_b,\\
u,\ \tau_x(u)n_x+\tau_z(u)n_z
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

代码默认反演例子取 $m=1$，此时底部项是

$$
\beta u\chi.
$$

若 `pde.f` 为空，`FirstOrderP2.m` 使用表面坡度给出驱动力

$$
f=-\rho g\,S_x.
$$

在 slab 中 $S_x=-\alpha$，因此当 $\rho,g>0$ 时 $f=\rho g\alpha$。

## 弱形式

设测试函数为 $\varphi$。正问题弱形式是求周期函数 $u$，使得

$$
R(u,\beta;\varphi)=0
\qquad\forall\varphi,
$$

其中

$$
R(u,\beta;\varphi)
=
\int_\Omega
\eta(u)
\left(
4u_x\varphi_x+u_z\varphi_z
\right)\,d\Omega
+
\int_{\Gamma_b}
\beta |u t_x|^{m-1}u\chi\,\varphi\,ds
-
\int_\Omega f\varphi\,d\Omega.
$$

对 $m=1$，

$$
R(u,\beta;\varphi)
=
\int_\Omega
\eta(u)
\left(
4u_x\varphi_x+u_z\varphi_z
\right)\,d\Omega
+
\int_{\Gamma_b}
\beta u\chi\,\varphi\,ds
-
\int_\Omega f\varphi\,d\Omega.
$$

`FirstOrderP2.m` 的 Picard 步固定当前速度中的黏度和底部系数，组装

$$
K_{\rm pic}(u^k)u^{k+1}=F.
$$

收敛后若 `assemble_tangent=true`，代码进一步组装一致切线
$K=\partial R/\partial u$，供伴随和 Gauss--Newton 使用。

## 目标函数

反问题使用顶部水平速度观测 $u_{\rm obs}$。连续目标函数可写为

$$
J(u)
=
\frac12
\int_{\Gamma_t}
\left(u-u_{\rm obs}\right)^2\,ds.
$$

`FOAdjInvSlabBed.m` 使用离散归一化版本：

$$
J_h(u)
=
\frac12
\frac{
\sum_{i\in\Gamma_t} w_i (u_i-u_{{\rm obs},i})^2
}{
\sum_{i\in\Gamma_t} w_i u_{{\rm obs},i}^2
}.
$$

其状态方向导数为

$$
\delta J(u)[\tilde u]
=
\int_{\Gamma_t}
(u-u_{\rm obs})\tilde u\,ds.
$$

离散代码中这对应

```matlab
observationGradient(topDof) = topWeight.*residual/dataNormSquared;
```

## 对数参数化

反演变量不是直接的 $\beta$，而是

$$
q=\log\beta,
\qquad
\beta=\exp(q).
$$

因此参数方向满足

$$
\delta\beta=\beta\,\delta q.
$$

`FOAdjInvSlabBed.m` 用周期 P1 函数表示 $q$ 和 $\beta$：

$$
\beta(x)=\operatorname{periodicP1}(x;\{x_j,\beta_j\}_{j=1}^{N_m}).
$$

## 线性化

令

$$
E(u,\varphi)=4u_x\varphi_x+u_z\varphi_z.
$$

黏度方向导数为

$$
\eta'(u)[\tilde u]
=
\eta(u)\,
\frac{1-n}{2n}\,
\frac{
2u_x\tilde u_x+\frac12u_z\tilde u_z
}{
\varepsilon_{\rm II}(u)+\varepsilon_{\rm reg}^2
}.
$$

因此正问题残差对状态的一阶变分为

$$
R_u(u,\beta)[\tilde u,\varphi]
=
\int_\Omega
\left[
\eta(u)E(\tilde u,\varphi)
+
\eta'(u)[\tilde u]E(u,\varphi)
\right]d\Omega
+
\int_{\Gamma_b}
\beta\chi\,\tilde u\,\varphi\,ds
$$

这里取 $m=1$。这就是 `assembleviscoustangent` 和 `assemblebedtangent` 的连续对应。

参数方向的变分为

$$
R_q(u,\beta)[\delta q,\varphi]
=
\int_{\Gamma_b}
\beta\,\delta q\,\chi\,u\,\varphi\,ds.
$$

这对应 `assemblebetadirection`；脚本把每个参数基函数方向组装成矩阵

$$
G_{ij}=R_q(u,\beta)[\psi_j,\phi_i].
$$

## 增量正问题

给定参数方向 $\delta q$，状态方向 $\tilde u$ 满足

$$
R_u(u,\beta)[\tilde u,\varphi]
=
-
R_q(u,\beta)[\delta q,\varphi]
\qquad\forall\varphi.
$$

在离散形式中，

$$
K\tilde U=-G\,\delta Q.
$$

这正是 Gauss--Newton 乘子中

```matlab
incrementalMaster = eqn.tangent\(-G*direction);
```

的含义。

## 伴随问题

伴随变量记为 $v$。为了消去目标函数方向导数中的 $\tilde u$，取 $v$ 满足

$$
R_u(u,\beta)[\varphi,v]
=
-
\delta J(u)[\varphi]
\qquad\forall\varphi.
$$

这里使用线性化算子的转置。有限维形式是

$$
K^T V=-C,
$$

其中 $C$ 是目标函数对状态的梯度。代码对应

```matlab
adjoint = eqn.tangent'\(-observationGradientMaster);
```

若把伴随方程写成强形式，可以理解为其顶部边界由当前观测误差驱动：

$$
\tau_x^*(v)n_x+\tau_z^*(v)n_z
=
-(u-u_{\rm obs})
\qquad\text{on }\Gamma_t.
$$

底部伴随边界项为

$$
\tau_x^*(v)n_x+\tau_z^*(v)n_z+\beta\chi v=0
\qquad\text{on }\Gamma_b
$$

这里取 $m=1$。其中 $\tau^*$ 表示一致切线算子的伴随作用。对这里的标量 FO
黏性切线，弱形式矩阵通常是对称的；代码仍显式使用 `eqn.tangent'`，保持与一般伴随
记号一致。

## 梯度公式

由增量正问题和伴随问题，

$$
\delta J(q)[\delta q]
=
R_q(u,\beta)[\delta q,v]
=
\int_{\Gamma_b}
\beta\,\delta q\,\chi\,u\,v\,ds.
$$

因此相对于 $q$ 的连续梯度密度是

$$
g_q=\beta\chi u v
\qquad\text{on }\Gamma_b.
$$

离散参数空间中，梯度为

$$
\nabla J_h(Q)=G^T V.
$$

代码对应

```matlab
gradient = G'*adjoint;
```

这个符号约定与 `adjoint = K'\(-C)` 配套使用。

## Gauss--Newton 方向

Gauss--Newton 方法求解

$$
\left(H_{\rm GN}+\lambda I\right)s=-\nabla J_h(Q),
$$

其中 $\lambda$ 是 Levenberg--Marquardt 阻尼。对给定方向 $d$，
`gaussnewtonproduct` 按下面三步计算
$H_{\rm GN}d+\lambda d$。

第一步，解增量正问题：

$$
K\tilde U=-Gd.
$$

第二步，把增量状态限制到顶部观测自由度，形成增量观测梯度

$$
\tilde C_i
=
\frac{w_i\tilde U_i}{\sum_j w_j u_{{\rm obs},j}^2},
\qquad i\in\Gamma_t.
$$

第三步，解增量伴随问题：

$$
K^T\tilde V=-\tilde C.
$$

最后投回参数空间：

$$
H_{\rm GN}d
=
G^T\tilde V.
$$

加上阻尼后得到

$$
\left(H_{\rm GN}+\lambda I\right)d
=
G^T\tilde V+\lambda d.
$$

代码对应

```matlab
incrementalMaster = eqn.tangent\(-G*direction);
incrementalAdjoint = eqn.tangent'\(-incrementalObservationMaster);
product = G'*incrementalAdjoint+lambda*direction;
```

## 逻辑关系

四个核心关系可以压缩为

$$
R(u,q)=0,
$$

$$
R_u(u,q)\tilde u=-R_q(u,q)\delta q,
$$

$$
R_u(u,q)^*v=-J_u(u),
$$

$$
H_{\rm GN}\delta q
=
R_q(u,q)^*\tilde v,
\qquad
R_u(u,q)^*\tilde v=-J_{uu}(u)\tilde u.
$$

在 `FOAdjInvSlabBed.m` 中：

- `solveforward` 解正问题；
- `eqn.tangent` 是 $R_u$ 的离散矩阵 $K$；
- `G` 是 $R_q$ 的离散矩阵；
- `adjoint = K'\(-C)` 解伴随问题；
- `gradient = G'*adjoint` 给出参数梯度；
- `gaussnewtonproduct` 实现 $H_{\rm GN}d+\lambda d$。

## 符号表

| 符号 | 含义 |
| :---: | :--- |
| $\Omega$ | 倾斜 slab 区域 |
| $\Gamma_t$ | 上表面 |
| $\Gamma_b$ | 底部 |
| $\Gamma_p$ | 左右周期边界 |
| $u$ | 水平速度，也是 FO 正问题主未知量 |
| $v$ | 伴随变量 |
| $\tilde u$ | 增量正问题状态方向 |
| $\tilde v$ | 增量伴随变量 |
| $A,n$ | Glen 流律参数 |
| $\eta(u)$ | 有效黏度 |
| $\varepsilon_{\rm II}(u)$ | FO 应变率不变量 |
| $\varepsilon_{\rm reg}$ | 黏度正则化参数 |
| $\beta$ | 底部摩擦系数 |
| $q=\log\beta$ | 反演参数 |
| $\chi=t_x^2$ | 底部切向几何因子 |
| $R$ | 正问题弱残差 |
| $R_u$ | 残差对状态的一阶导数，也就是一致切线 |
| $R_q$ | 残差对对数参数的一阶导数 |
| $G$ | 离散参数导数矩阵 |
| $J$ | 顶部速度失配目标函数 |
| $K$ | `eqn.tangent`，离散一致切线矩阵 |
| $\lambda$ | Levenberg--Marquardt 阻尼参数 |
