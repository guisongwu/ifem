# 非线性 Stokes 反问题中的四个方程

本文从连续 PDE 角度整理非线性 Stokes 冰流反问题中的四类方程：

1. 正问题；
2. 增量正问题；
3. 伴随问题；
4. 增量伴随问题。

目标是说明这些方程从哪里来、每个符号表示什么，以及它们如何组成伴随反演和 Gauss--Newton 方法。本文只讨论连续形式，不涉及离散矩阵或数值实现。

本文取底部线性滑移，即滑移指数 $m=1$。通用 $m$ 的版本见 `stress-linearization.md`。

## 变分与方向导数

本文中的一阶变分可以理解为 Gateaux 方向导数。设 $F$ 是定义在某个函数空间上的泛函或算子，$u$ 是当前点，$\tilde u$ 是扰动方向，则

$$
\delta F(u)[\tilde u]
=
\left.\frac{\mathrm{d}}{\mathrm{d}t}F(u+t\tilde u)\right|_{t=0}
=
\lim_{t\to 0}
\frac{F(u+t\tilde u)-F(u)}{t}.
$$

这里 $\delta F(u)$ 表示 $F$ 在 $u$ 处的一阶线性化，方括号中的 $\tilde u$ 表示这个线性化算子作用的方向。

若 $F$ 是 Frechet 可微的，则还可以写成

$$
\delta F(u)[\tilde u]=F'(u)\tilde u,
$$

其中 $F'(u)$ 是 $F$ 在 $u$ 处的导数算子。

本文后面出现的目标函数变分、残差变分、应力变分都采用以上记号。

## 记号约定

无修饰变量表示当前正问题或伴随问题的主变量：

- $(\boldsymbol u,p)$ 表示给定 $\beta$ 后的正问题速度和压力；
- $(\boldsymbol v,r)$ 表示由当前顶部观测误差驱动的伴随速度和伴随压力。

带 $\tilde{\ }$ 的变量表示增量问题中的未知量：

- $(\tilde{\boldsymbol u},\tilde p)$ 表示参数扰动引起的状态增量；
- $(\tilde{\boldsymbol v},\tilde r)$ 表示对应的增量伴随变量。

$\delta$ 表示变分、方向导数或参数扰动方向。如

- 若函数 $F$ 对变量 $x$ 求变分并作用在方向 $h$ 上，本文记作 $\delta F(x)[h]$, 有时也记作 $\delta_xF[h]$；
- 参数扰动仍记作 $\delta\beta,\delta q$;
- $\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]$ 表示应力 $\boldsymbol\sigma$ 对状态变量 $(\boldsymbol u,p)$ 的一阶变分，作用在方向 $(\tilde{\boldsymbol u},\tilde p)$ 上。

## 几何、应力与滑移律

设冰体区域为 $\Omega\subset\mathbb R^2$，边界分为

$$
\partial\Omega
=
\Gamma_t\cup\Gamma_b\cup\Gamma_p.
$$

其中 $\Gamma_t$ 是上表面，$\Gamma_b$ 是底部，$\Gamma_p$ 表示左右周期边界。记 $\boldsymbol n$ 为区域外法向，$\boldsymbol e_x=(1,0)^T$ 为水平方向单位向量。

底部切向投影算子定义为

$$
\boldsymbol T
=
\boldsymbol I-\boldsymbol n\otimes\boldsymbol n.
$$

这里的张量积和切向投影记号已在 `tensor-notation.md` 中说明。任意速度 $\boldsymbol w$ 的底部切向分量记为

$$
\boldsymbol w_t=\boldsymbol T\boldsymbol w.
$$

速度 $\boldsymbol w$ 的应变率张量记为 $\dot{\boldsymbol\varepsilon}_{\boldsymbol w}$，定义为

$$
\dot{\boldsymbol\varepsilon}_{\boldsymbol w}
=
\frac12
\left(
\nabla\boldsymbol w+\nabla\boldsymbol w^T
\right).
$$

非线性 Stokes 应力为

$$
\boldsymbol\sigma
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
-
p\boldsymbol I.
$$

这里 $\eta(\boldsymbol u)$ 是有效黏度。Glen 型流律中，$\eta$ 依赖应变率不变量，因此间接依赖速度 $\boldsymbol u$。这也是正问题非线性的主要来源。

本文取 $m=1$ 的线性底部滑移律，底部摩擦力为

$$
\beta\boldsymbol u_t.
$$

于是底部边界上同时有不可穿透条件和切向滑移条件：

$$
\boldsymbol u\cdot\boldsymbol n=0,
\qquad
\boldsymbol T\boldsymbol\sigma\boldsymbol n
+
\beta\boldsymbol u_t
=
\boldsymbol 0
\qquad\text{on }\Gamma_b.
$$

## 正问题

给定底部摩擦系数 $\beta$，正问题是求 $(\boldsymbol u,p)$，满足

$$
\left\{
\begin{aligned}
-\nabla\cdot\boldsymbol\sigma &= \textcolor{red}{\rho\boldsymbol g}
&&\text{in }\Omega,\\
\nabla\cdot\boldsymbol u &=0
&&\text{in }\Omega,\\
\boldsymbol\sigma\boldsymbol n &= \boldsymbol 0
&&\text{on }\Gamma_t,\\
\boldsymbol u\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T\boldsymbol\sigma\boldsymbol n
+
\beta\boldsymbol u_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol u,\ \boldsymbol\sigma\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

其中

$$
\boldsymbol\sigma
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
-
p\boldsymbol I.
$$

其中 $\boldsymbol\sigma_{\rm visc}=2\eta(\boldsymbol u)\dot{\boldsymbol\varepsilon}_{\boldsymbol u}$ 是黏性应力，也可以把正问题抽象写成非线性残差方程:

$$
\mathcal R(\boldsymbol u,p,\beta)=0.
$$

这里 $\mathcal R$ 收集了动量方程、不可压缩条件、顶部自然边界条件、底部不可穿透条件、底部滑移条件和周期条件的残差。若 $(\boldsymbol u,p)$ 是给定 $\beta$ 下的正问题解，则这些残差同时为零。

## 目标函数

假设顶部有水平速度观测 $u_{\rm obs}$, 目标函数取为

$$
J(\boldsymbol u)
=
\frac12
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)^2\,\mathrm{d}s,
$$

其中

$$
u_x=\boldsymbol u\cdot\boldsymbol e_x.
$$

对状态变量 $\boldsymbol u$ 沿方向 $\tilde{\boldsymbol u}$ 求一阶变分。

$$
\begin{aligned}
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
&=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)\tilde u_x\,ds \\
&=
\int_{\Gamma_t}
\left[
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\right]
\cdot
\tilde{\boldsymbol u}\,\mathrm{d}s.
\end{aligned}
$$

这个顶部边界项会成为伴随方程顶部边界条件的来源。

## 线性化准备

反问题需要知道参数变化如何影响目标函数。设参数有扰动 $\delta\beta$，由正问题引起的状态扰动记为 $(\tilde{\boldsymbol u},\tilde p)$。对残差方程

$$
\mathcal R(\boldsymbol u,p,\beta)=0
$$

做一阶变分，得到

$$
\delta\mathcal R(\boldsymbol u,p,\beta)
[\tilde{\boldsymbol u},\tilde p,0]
+
\delta\mathcal R(\boldsymbol u,p,\beta)
[0,0,\delta\beta]
=
0.
$$

其中第一项是残差对状态变量 $(\boldsymbol u,p)$ 的一阶变分，第二项是残差对参数 $\beta$ 的一阶变分。

### 应力线性化

考虑扰动状态

$$
\boldsymbol u_\epsilon
=
\boldsymbol u+\epsilon\tilde{\boldsymbol u},
\qquad
p_\epsilon
=
p+\epsilon\tilde p.
$$

应力的一阶变化为

$$
\boxed{
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
+
2\eta'(\boldsymbol u)[\tilde{\boldsymbol u}]
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
-
\tilde p\boldsymbol I
}
$$

这三个项分别来自应变率变化、黏度变化和压力变化。其中

$$
\eta(\bm u)
=
\frac12 A^{-1/n}\varepsilon_\text{II}^{\frac{1-n}{2n}}.
$$

也可以把黏性部分写成四阶切线张量形式。

若黏度只依赖应变率第二不变量

$$
\varepsilon_\text{II}
=
\frac12
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\boldsymbol u},
$$

则

$$
\eta'(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\frac{\partial\eta}{\partial\varepsilon_\text{II}}
\left(
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right).
$$

因此黏性应力的一阶变化可以写成

$$
\delta\boldsymbol\sigma_{\rm visc}
=
2\eta
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
+
2
\frac{\partial\eta}{\partial\varepsilon_\text{II}}
\left(
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}.
$$

在当前速度场 $\boldsymbol u$ 固定时，上式是关于 $\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$ 的线性映射，因此可以用四阶切线张量 $\mathbb C(\boldsymbol u)$ 表示。相关张量记号和分量形式已在 `tensor-notation.md` 中说明。这里直接记为

$$
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
=
2\eta
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
+
2
\frac{\partial\eta}{\partial\varepsilon_\text{II}}
\left(
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}.
$$

于是完整线性化应力可以写成

$$
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
-
\tilde p\boldsymbol I.
$$

若只保留 $2\eta\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$，就是冻结黏度的 Picard 线性化；包含黏度导数项的完整表达式才是一致切线。

### 底部滑移线性化

因为本文取 $m=1$，底部摩擦项是

$$
\beta\boldsymbol u_t.
$$

对状态方向 $\tilde{\boldsymbol u}$，有

$$
\delta_{\boldsymbol u}
\left(
\beta\boldsymbol u_t
\right)
[\tilde{\boldsymbol u}]
=
\beta\tilde{\boldsymbol u}_t.
$$

对参数方向 $\delta\beta$，有

$$
\delta_\beta
\left(
\beta\boldsymbol u_t
\right)
[\delta\beta]
=
\delta\beta\boldsymbol u_t.
$$

如果使用对数参数

$$
q=\log\beta,
$$

则

$$
\delta\beta=\beta\,\delta q.
$$

## 增量正问题

增量正问题描述：参数沿方向 $\delta\beta$ 改变时，状态的一阶响应 $(\tilde{\boldsymbol u},\tilde p)$ 满足什么方程。

由残差线性化

$$
\delta\mathcal R(\boldsymbol u,p,\beta)
[\tilde{\boldsymbol u},\tilde p,0]
=
-
\delta\mathcal R(\boldsymbol u,p,\beta)
[0,0,\delta\beta],
$$

得到增量正问题的强形式：

$$
\left\{
\begin{aligned}
-\nabla\cdot
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\nabla\cdot\tilde{\boldsymbol u} &=0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
&=\boldsymbol 0
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol u}\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
+
\beta\tilde{\boldsymbol u}_t
&=
\textcolor{red}{-\delta\beta\boldsymbol u_t}
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol u},\
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

源项只出现在底部，是因为 $\beta$ 只出现在底部滑移边界条件中。

## 伴随问题

如果直接计算约化目标函数在方向 $\delta\beta$ 上的一阶变分 $\delta J(\beta)[\delta\beta]$，需要知道每个参数方向 $\delta\beta$ 对应的状态扰动 $\tilde{\boldsymbol u}$。伴随方法的目的，是用一个伴随问题消去 $\tilde{\boldsymbol u}$，从而直接得到关于 $\delta\beta$ 的梯度表达。

伴随问题来自增量正问题线性化算子的转置，以及顶部目标函数变分与边界项的抵消。完整的分部积分推导见 [adjoint.md](adjoint.md)。本文只记录所采用的符号约定和最终强形式。

记增量正问题的状态线性化算子为

$$
\mathcal A(\tilde{\boldsymbol u},\tilde p)
=
\delta\mathcal R(\boldsymbol u,p,\beta)
[\tilde{\boldsymbol u},\tilde p,0].
$$

伴随算子 $\mathcal A^*$ 由同一弱形式配对定义：

$$
\left\langle
\mathcal A(\tilde{\boldsymbol u},\tilde p),
(\boldsymbol v,r)
\right\rangle
=
\left\langle
(\tilde{\boldsymbol u},\tilde p),
\mathcal A^*(\boldsymbol v,r)
\right\rangle
+
\text{boundary terms}.
$$

若增量正问题的线性化应力为

$$
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
-
\tilde p\boldsymbol I,
$$

则伴随线性化应力定义为

$$
\delta\boldsymbol\sigma^*(\bm u,p)[\boldsymbol v,r]
=
\mathbb C(\boldsymbol u)^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
-
r\boldsymbol I.
$$

这里 $\mathbb C(\boldsymbol u)^T$ 是相对于张量内积的转置，具体定义和指标交换关系已在 `tensor-notation.md` 中说明。对本文中的 Glen 型黏性切线，有

$$
\mathbb C(\boldsymbol u)^T=\mathbb C(\boldsymbol u).
$$

这时体内的伴随线性化应力在形式上和增量正问题中的线性化应力相同，只是方向变量换成 $(\boldsymbol v,r)$。

伴随问题为求 $(\boldsymbol v,r)$，满足

$$
\left\{
\begin{aligned}
-\nabla\cdot
\delta\boldsymbol\sigma^*(\bm u,p)[\boldsymbol v,r]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\nabla\cdot\boldsymbol v &=0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma^*(\bm u,p)[\boldsymbol v,r]\boldsymbol n
&=
\textcolor{red}{-\left(u_x-u_{\rm obs}\right)\boldsymbol e_x}
&&\text{on }\Gamma_t,\\
\boldsymbol v\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma^*(\bm u,p)[\boldsymbol v,r]\boldsymbol n
+
\beta\boldsymbol v_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol v,\
\delta\boldsymbol\sigma^*(\bm u,p)[\boldsymbol v,r]\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

底部项中出现 $\beta\boldsymbol v_t$，是因为 $m=1$ 时底部滑移线性化算子为 $\beta\boldsymbol T$，而 $\boldsymbol T$ 是自伴随投影。

## 梯度公式

由增量正问题和伴随问题相消，可以把目标函数对参数方向的变化写成底部边界积分：

$$
\delta J(\beta)[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

如果使用 $q=\log\beta$，则 $\delta\beta=\beta\delta q$，于是

$$
\delta J(q)[\delta q]
=
\int_{\Gamma_b}
\beta\delta q\,
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

因此相对于 $q$ 的连续梯度密度为

$$
g_q
=
\beta
\boldsymbol u_t\cdot\boldsymbol v_t.
$$

符号可能随伴随变量的符号约定整体改变；只要伴随方程和梯度公式使用同一约定，下降方向是一致的。

## 增量伴随问题

增量伴随问题用于计算梯度沿参数方向 $\delta\beta$ 的变化，也就是 Hessian-vector product。

给定参数方向 $\delta\beta$，先解增量正问题得到状态方向

$$
(\tilde{\boldsymbol u},\tilde p).
$$

目标函数顶部源项

$$
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
$$

沿状态方向的一阶变化是

$$
\tilde u_x\boldsymbol e_x.
$$

因此增量伴随问题的顶部边界条件由 $\tilde u_x$ 驱动。Gauss--Newton 近似忽略正问题二阶导数和伴随方程中由非线性算子二阶导数带来的项，只保留观测算子的二阶项。于是增量伴随问题为求 $(\tilde{\boldsymbol v},\tilde r)$，满足

$$
\left\{
\begin{aligned}
-\nabla\cdot
\delta\boldsymbol\sigma^*(\bm u,p)[\tilde{\boldsymbol v},\tilde r]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\nabla\cdot\tilde{\boldsymbol v} &=0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma^*(\bm u,p)[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
&=
\textcolor{red}{-\tilde u_x\boldsymbol e_x}
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol v}\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma^*(\bm u,p)[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
+
\beta\tilde{\boldsymbol v}_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol v},\
\delta\boldsymbol\sigma^*(\bm u,p)[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

它和伴随问题的区别只在顶部边界源项：

$$
-\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\quad\longrightarrow\quad
-\tilde u_x\boldsymbol e_x.
$$

前者由当前观测误差驱动，后者由增量状态在观测量上的变化驱动。

## Gauss--Newton Hessian-vector product

Gauss--Newton 方法需要计算 $H_{\rm GN}\delta q$。连续层面可以按下面的顺序理解。

首先给定参数方向 $\delta q$，令

$$
\delta\beta=\beta\delta q.
$$

然后解增量正问题，得到 $(\tilde{\boldsymbol u},\tilde p)$。这个解给出顶部观测量的一阶变化 $\tilde u_x$。

接着解增量伴随问题，得到 $(\tilde{\boldsymbol v},\tilde r)$。最后把信息投回到底部参数空间：对任意测试方向 $\zeta$，

$$
\left\langle
H_{\rm GN}\delta q,\zeta
\right\rangle
=
\int_{\Gamma_b}
\beta
\left(
\boldsymbol u_t\cdot\tilde{\boldsymbol v}_t
\right)
\zeta\,ds.
$$

也就是说，$H_{\rm GN}\delta q$ 是底部函数

$$
\beta
\boldsymbol u_t\cdot\tilde{\boldsymbol v}_t
$$

在参数空间内积下的表示或投影。这个表达式保留了观测残差对参数方向的二阶影响，但忽略了完整 Hessian 中含有正问题二阶导数的部分，所以它是 Gauss--Newton 近似。

## 四个方程的逻辑关系

四个方程可以压缩成下面的抽象关系。

正问题：

$$
\mathcal R(\boldsymbol u,p,\beta)=0.
$$

增量正问题：

$$
\delta\mathcal R(\boldsymbol u,p,\beta)
[\tilde{\boldsymbol u},\tilde p,0]
=
-
\delta\mathcal R(\boldsymbol u,p,\beta)
[0,0,\delta\beta].
$$

伴随问题：

$$
\delta\mathcal R(\boldsymbol u,p,\beta)^*
[\boldsymbol v,r]
=
-
\delta J(\boldsymbol u)[\cdot].
$$

增量伴随问题：

$$
\delta\mathcal R(\boldsymbol u,p,\beta)^*
[\tilde{\boldsymbol v},\tilde r]
=
-
\delta^2J(\boldsymbol u)
[\tilde{\boldsymbol u},\cdot].
$$

这里 $\delta J(\boldsymbol u)[\cdot]$ 是目标函数对状态的一阶变分，$\delta^2J(\boldsymbol u)[\tilde{\boldsymbol u},\cdot]$ 是其二阶变分在第一方向固定为 $\tilde{\boldsymbol u}$ 后得到的线性泛函。

对于顶部水平速度观测，有

$$
\delta J(\boldsymbol u)[\boldsymbol w]
=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)w_x\,ds,
$$

以及

$$
\delta^2J(\boldsymbol u)
[\tilde{\boldsymbol u},\boldsymbol w]
=
\int_{\Gamma_t}
\tilde u_x w_x\,ds.
$$

## 伴随反演思想小结

伴随反演的核心是避免对每一个参数自由度都单独求一次状态导数。正问题给出映射

$$
\beta \longmapsto \boldsymbol u(\beta).
$$

反问题希望通过顶部观测误差来调整底部参数 $\beta$。如果直接计算梯度，需要知道每个参数扰动 $\delta\beta$ 引起的速度扰动 $\tilde{\boldsymbol u}$，这会导致大量增量正问题求解。

伴随方法的做法是先把正问题在当前解处线性化：

$$
\mathcal A(\tilde{\boldsymbol u},\tilde p)
=
-
\delta\mathcal R(\boldsymbol u,p,\beta)
[0,0,\delta\beta].
$$

然后引入伴随变量 $(\boldsymbol v,r)$，令它满足线性化算子的转置问题：

$$
\mathcal A^*(\boldsymbol v,r)
=
-
\delta J(\boldsymbol u)[\cdot].
$$

这样目标函数变分中的状态扰动 $\tilde{\boldsymbol u}$ 可以被消去，最终梯度只表现为底部边界上的参数项：

$$
\delta J(\beta)[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

因此，一次梯度计算只需要一次正问题和一次伴随问题。进一步地，Gauss--Newton 方法需要 Hessian-vector product；这时对给定方向先解增量正问题，再解增量伴随问题，就可以得到该方向上的二阶信息。

简而言之：

- 正问题把底部参数 $\beta$ 传到顶部速度；
- 伴随问题把顶部观测误差传回到底部；
- 梯度由底部速度、伴随速度和滑移律共同给出；
- 增量正问题和增量伴随问题用于计算 Gauss--Newton Hessian 对方向的作用。

## 符号表

| 符号 | 含义 |
| :---: | :---: |
| $\Omega$ | 冰体区域 |
| $\Gamma_t$ | 上表面 |
| $\Gamma_b$ | 底部 |
| $\Gamma_p$ | 周期边界 |
| $\boldsymbol n$ | 外法向 |
| $\boldsymbol T=\boldsymbol I-\boldsymbol n\otimes\boldsymbol n$ | 切向投影 |
| $\boldsymbol u$ | 正问题速度 |
| $p$ | 正问题压力 |
| $\boldsymbol v$ | 伴随速度 |
| $r$ | 伴随压力 |
| $\tilde{\boldsymbol u}$ | 增量正问题速度方向 |
| $\tilde p$ | 增量正问题压力方向 |
| $\tilde{\boldsymbol v}$ | 增量伴随速度 |
| $\tilde r$ | 增量伴随压力 |
| $\beta$ | 底部摩擦系数 |
| $q=\log\beta$ | 对数参数 |
| $\delta\beta$ | 摩擦系数扰动 |
| $\delta q$ | 对数参数扰动 |
| $\dot{\boldsymbol\varepsilon}_{\boldsymbol u}$ | 正问题速度的应变率张量 |
| $\eta(\boldsymbol u)$ | 有效黏度 |
| $\boldsymbol\sigma$ | 正问题应力 |
| $\delta\boldsymbol\sigma(\bm u,p)[\cdot,\cdot]$ | 应力在 $(\bm u,p)$ 处的一阶线性化 |
| $\delta\boldsymbol\sigma^*(\bm u,p)[\cdot,\cdot]$ | 线性化应力算子的伴随 |
| $u_{\rm obs}$ | 顶部水平速度观测 |
| $J$ | 反问题目标函数 |
| $\mathbb C(\boldsymbol u)$ | 黏性一致切线张量 |
