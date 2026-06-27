# 非线性 Stokes 反问题中的四个方程

本文从连续数学角度整理非线性 Stokes 冰流反问题中的四类方程：

1. 正问题；
2. 增量正问题；
3. 伴随问题；
4. 增量伴随问题。

目标是说明这些方程各自从哪里来、每个符号表示什么，以及它们如何组成伴随反演和 Gauss--Newton 方法。本文只讨论连续形式，不涉及离散矩阵或数值实现。

本文取底部线性滑移，即滑移指数 $m=1$。通用 $m$ 的版本见 `stress-linearization.md`。

## 记号约定

无修饰变量表示当前正问题或伴随问题的主变量：

- $(\boldsymbol u,p)$ 表示给定 $\beta$ 后的正问题速度和压力；
- $(\boldsymbol v,r)$ 表示由当前顶部观测误差驱动的伴随速度和伴随压力。

带 $\tilde{\ }$ 的变量表示增量问题中的未知量：

- $(\tilde{\boldsymbol u},\tilde p)$ 表示参数扰动引起的状态增量；
- $(\tilde{\boldsymbol v},\tilde r)$ 表示对应的增量伴随变量。

$\delta$ 表示变分、方向导数或参数扰动方向，例如

$$
\delta\beta,\qquad \delta q,\qquad \delta J.
$$

特别地，

$$
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]
$$

表示应力 $\boldsymbol\sigma$ 对状态变量 $(\boldsymbol u,p)$ 的一阶变分，作用在方向 $(\tilde{\boldsymbol u},\tilde p)$ 上。这里的 $\delta\boldsymbol\sigma$ 不是新的物理应力，而是线性化应力算子。

## 几何与边界

设冰体区域为 $\Omega\subset\mathbb R^2$，边界分为

$$
\partial\Omega
=
\Gamma_t\cup\Gamma_b\cup\Gamma_p.
$$

其中：

- $\Gamma_t$ 是上表面；
- $\Gamma_b$ 是底部；
- $\Gamma_p$ 表示左右周期边界；
- $\boldsymbol n$ 是区域外法向；
- $\boldsymbol e_x=(1,0)^T$ 是水平方向单位向量。

底部切向投影算子定义为

$$
\boldsymbol T
=
\boldsymbol I-\boldsymbol n\otimes\boldsymbol n.
$$

因此任意速度 $\boldsymbol u$ 的底部切向分量为

$$
\boldsymbol u_t
=
\boldsymbol T\boldsymbol u.
$$

同理，

$$
\boldsymbol v_t=\boldsymbol T\boldsymbol v,\qquad
\tilde{\boldsymbol u}_t=\boldsymbol T\tilde{\boldsymbol u},\qquad
\tilde{\boldsymbol v}_t=\boldsymbol T\tilde{\boldsymbol v}.
$$

底部不可穿透条件是

$$
\boldsymbol u\cdot\boldsymbol n=0
\qquad\text{on }\Gamma_b.
$$

## 应力与滑移律

速度应变率张量为

$$
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
=
\frac12
\left(
\nabla\boldsymbol u+\nabla\boldsymbol u^T
\right).
$$

非线性 Stokes 应力为

$$
\boldsymbol\sigma
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
-
p\boldsymbol I.
$$

这里 $\eta(\boldsymbol u)$ 是有效黏度。Glen 型流律中，$\eta$ 依赖应变率不变量，因此间接依赖速度 $\boldsymbol u$。这也是正问题非线性的主要来源。

本文取 $m=1$ 的线性底部滑移律，底部摩擦力为

$$
\beta\boldsymbol u_t.
$$

于是底部滑移边界条件为

$$
\boldsymbol T\boldsymbol\sigma\boldsymbol n
+
\beta\boldsymbol u_t
=
\boldsymbol 0
\qquad\text{on }\Gamma_b.
$$

第一项是底部牵引的切向分量，第二项是底部摩擦力。

## 正问题

给定底部摩擦系数 $\beta$，正问题是求 $(\boldsymbol u,p)$，满足

$$
\left\{
\begin{aligned}
\nabla\cdot\boldsymbol u &=0
&&\text{in }\Omega,\\
-\nabla\cdot\boldsymbol\sigma &= \rho\boldsymbol g
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
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
-
p\boldsymbol I.
$$

也可以把正问题抽象写成非线性残差方程

$$
\mathcal R(\boldsymbol u,p,\beta)=0.
$$

这里 $\mathcal R$ 收集了动量方程、不可压缩条件、顶部自然边界条件、底部不可穿透条件、底部滑移条件和周期条件的残差。

例如体内动量残差是

$$
-\nabla\cdot\boldsymbol\sigma-\rho\boldsymbol g,
$$

底部滑移残差是

$$
\boldsymbol T\boldsymbol\sigma\boldsymbol n+\beta\boldsymbol u_t.
$$

若 $(\boldsymbol u,p)$ 是给定 $\beta$ 下的正问题解，则这些残差同时为零。

## 目标函数

假设顶部有水平速度观测 $u_{\rm obs}$。目标函数取为

$$
J(\boldsymbol u)
=
\frac12
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)^2\,ds,
$$

其中

$$
u_x=\boldsymbol u\cdot\boldsymbol e_x.
$$

目标函数的一阶变分为

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)\tilde u_x\,ds.
$$

这里出现 $\tilde u_x$ 的原因是：一阶变分是在状态变量
$\boldsymbol u$ 的扰动方向 $\tilde{\boldsymbol u}$ 上求方向导数。
令

$$
\boldsymbol u_\epsilon
=
\boldsymbol u+\epsilon\tilde{\boldsymbol u},
$$

则水平速度分量随之变为

$$
(u_\epsilon)_x
=
u_x+\epsilon\tilde u_x,
$$

其中

$$
\tilde u_x=\tilde{\boldsymbol u}\cdot\boldsymbol e_x.
$$

因此

$$
J(\boldsymbol u+\epsilon\tilde{\boldsymbol u})
=
\frac12
\int_{\Gamma_t}
\left(u_x+\epsilon\tilde u_x-u_{\rm obs}\right)^2\,ds.
$$

对 $\epsilon$ 求导并令 $\epsilon=0$，得到

$$
\left.
\frac{d}{d\epsilon}
J(\boldsymbol u+\epsilon\tilde{\boldsymbol u})
\right|_{\epsilon=0}
=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)\tilde u_x\,ds.
$$

所以 $\tilde u_x$ 不是新的未知量，而是状态扰动
$\tilde{\boldsymbol u}$ 的水平分量。

由于

$$
\tilde u_x=\tilde{\boldsymbol u}\cdot\boldsymbol e_x,
$$

也可以写成

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
\left[
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\right]
\cdot
\tilde{\boldsymbol u}\,ds.
$$

这个顶部边界项会成为伴随方程顶部边界条件的来源。

## 一阶线性化

反问题需要知道参数变化如何影响目标函数。设参数有扰动 $\delta\beta$，由正问题引起的状态扰动记为 $(\tilde{\boldsymbol u},\tilde p)$。对残差方程

$$
\mathcal R(\boldsymbol u,p,\beta)=0
$$

做一阶变分，得到

$$
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p]
+
\mathcal R_\beta[\delta\beta]
=
0.
$$

其中 $\mathcal R_{(\boldsymbol u,p)}$ 是正问题残差对状态变量的一阶导数，$\mathcal R_\beta$ 是对参数的一阶导数。

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

应力的一阶变化定义为

$$
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]
=
\left.
\frac{d}{d\epsilon}
\boldsymbol\sigma(\boldsymbol u_\epsilon,p_\epsilon)
\right|_{\epsilon=0}.
$$

逐项求导得到

$$
\boxed{
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u})
+
2\eta'(\boldsymbol u)[\tilde{\boldsymbol u}]
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
-
\tilde p\boldsymbol I
}
$$

其中

$$
\eta'(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\left.
\frac{d}{d\epsilon}
\eta(\boldsymbol u+\epsilon\tilde{\boldsymbol u})
\right|_{\epsilon=0}.
$$

这三个项分别来自应变率变化、黏度变化和压力变化。

<!-- 如果把黏性部分写成四阶切线张量形式，则 -->

<!-- $$ -->
<!-- \delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p] -->
<!-- = -->
<!-- \mathbb C(\boldsymbol u): -->
<!-- \dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u}) -->
<!-- - -->
<!-- \tilde p\boldsymbol I. -->
<!-- $$ -->

<!-- 这里 $\mathbb C(\boldsymbol u)$ 是黏性一致切线张量。 -->

<!-- 下面说明这个四阶张量形式是怎么来的。记当前应变率和扰动应变率为 -->

<!-- $$ -->
<!-- \boldsymbol E -->
<!-- = -->
<!-- \dot{\boldsymbol\varepsilon}(\boldsymbol u), -->
<!-- \qquad -->
<!-- \tilde{\boldsymbol E} -->
<!-- = -->
<!-- \dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u}). -->
<!-- $$ -->

<!-- 黏性应力部分是 -->

<!-- $$ -->
<!-- \boldsymbol\sigma_{\rm visc} -->
<!-- = -->
<!-- 2\eta(\boldsymbol E)\boldsymbol E. -->
<!-- $$ -->

<!-- 给速度一个扰动 -->

<!-- $$ -->
<!-- \boldsymbol u_\epsilon -->
<!-- = -->
<!-- \boldsymbol u+\epsilon\tilde{\boldsymbol u}, -->
<!-- $$ -->

<!-- 则应变率也相应扰动为 -->

<!-- $$ -->
<!-- \boldsymbol E_\epsilon -->
<!-- = -->
<!-- \boldsymbol E+\epsilon\tilde{\boldsymbol E}. -->
<!-- $$ -->

<!-- 因此 -->

<!-- $$ -->
<!-- \boldsymbol\sigma_{\rm visc}(\boldsymbol u_\epsilon) -->
<!-- = -->
<!-- 2\eta(\boldsymbol E_\epsilon)\boldsymbol E_\epsilon. -->
<!-- $$ -->

<!-- 对 $\epsilon$ 求导并令 $\epsilon=0$，得到 -->

<!-- $$ -->
<!-- \left. -->
<!-- \frac{d}{d\epsilon} -->
<!-- \left( -->
<!-- 2\eta(\boldsymbol E_\epsilon)\boldsymbol E_\epsilon -->
<!-- \right) -->
<!-- \right|_{\epsilon=0} -->
<!-- = -->
<!-- 2\eta(\boldsymbol E)\tilde{\boldsymbol E} -->
<!-- + -->
<!-- 2\eta'(\boldsymbol E)[\tilde{\boldsymbol E}] -->
<!-- \boldsymbol E. -->
<!-- $$ -->

<!-- 第一项来自 $\boldsymbol E$ 本身的变化；第二项来自黏度 -->
<!-- $\eta$ 随应变率变化而产生的变化。 -->

<!-- 若黏度只依赖应变率第二不变量 -->

<!-- $$ -->
<!-- \varepsilon_{II} -->
<!-- = -->
<!-- \frac12\boldsymbol E:\boldsymbol E, -->
<!-- $$ -->

<!-- 则 -->

<!-- $$ -->
<!-- \eta'(\boldsymbol E)[\tilde{\boldsymbol E}] -->
<!-- = -->
<!-- \frac{\partial\eta}{\partial\varepsilon_{II}} -->
<!-- \left( -->
<!-- \boldsymbol E:\tilde{\boldsymbol E} -->
<!-- \right). -->
<!-- $$ -->

<!-- 于是黏性应力的一阶变化可以写成 -->

<!-- $$ -->
<!-- \delta\boldsymbol\sigma_{\rm visc} -->
<!-- = -->
<!-- 2\eta\tilde{\boldsymbol E} -->
<!-- + -->
<!-- 2 -->
<!-- \frac{\partial\eta}{\partial\varepsilon_{II}} -->
<!-- \left( -->
<!-- \boldsymbol E:\tilde{\boldsymbol E} -->
<!-- \right) -->
<!-- \boldsymbol E. -->
<!-- $$ -->

<!-- 这个表达式对扰动应变率 $\tilde{\boldsymbol E}$ 是线性的，因此可以定义四阶张量 -->
<!-- $\mathbb C(\boldsymbol u)$，使得 -->

<!-- $$ -->
<!-- \mathbb C(\boldsymbol u):\tilde{\boldsymbol E} -->
<!-- = -->
<!-- 2\eta\tilde{\boldsymbol E} -->
<!-- + -->
<!-- 2 -->
<!-- \frac{\partial\eta}{\partial\varepsilon_{II}} -->
<!-- \left( -->
<!-- \boldsymbol E:\tilde{\boldsymbol E} -->
<!-- \right) -->
<!-- \boldsymbol E. -->
<!-- $$ -->

<!-- 用分量形式写，就是 -->

<!-- $$ -->
<!-- C_{ijkl}(\boldsymbol u) -->
<!-- = -->
<!-- 2\eta\,\delta_{ik}\delta_{jl} -->
<!-- + -->
<!-- 2 -->
<!-- \frac{\partial\eta}{\partial\varepsilon_{II}} -->
<!-- E_{ij}E_{kl}, -->
<!-- $$ -->

<!-- 并且 -->

<!-- $$ -->
<!-- \left( -->
<!-- \mathbb C(\boldsymbol u):\tilde{\boldsymbol E} -->
<!-- \right)_{ij} -->
<!-- = -->
<!-- C_{ijkl}(\boldsymbol u)\tilde E_{kl}. -->
<!-- $$ -->

<!-- 由于 $\tilde{\boldsymbol E}$ 是对称应变率张量，这里的四阶张量只需要理解为作用在对称二阶张量空间上的线性算子。写成 $\mathbb C:\tilde{\boldsymbol E}$ 的好处是把“冻结黏度项”和“黏度导数项”合并成一个统一的线性化黏性算子： -->

<!-- $$ -->
<!-- \delta\boldsymbol\sigma_{\rm visc} -->
<!-- = -->
<!-- \mathbb C(\boldsymbol u):\tilde{\boldsymbol E}. -->
<!-- $$ -->

<!-- 再加上压力项的一阶变化 $-\tilde p\boldsymbol I$，就得到 -->

<!-- $$ -->
<!-- \delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p] -->
<!-- = -->
<!-- \mathbb C(\boldsymbol u): -->
<!-- \dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u}) -->
<!-- - -->
<!-- \tilde p\boldsymbol I. -->
<!-- $$ -->

<!-- 因此 $\mathbb C(\boldsymbol u)$ 不是额外引入的物理模型，而是非线性黏性应力 -->
<!-- $2\eta(\boldsymbol u)\dot{\boldsymbol\varepsilon}(\boldsymbol u)$ -->
<!-- 对当前速度场的一阶导数。若只保留 -->
<!-- $2\eta\tilde{\boldsymbol E}$，就是冻结黏度的 Picard 线性化；包含 -->
<!-- $\eta'(\boldsymbol E)[\tilde{\boldsymbol E}]$ 的完整表达式才是一致切线。 -->

### $\mathbb C(\boldsymbol u)$ 的显式形式与代码装配

记当前应变率和增量应变率为

$$
\boldsymbol E
=
\dot{\boldsymbol\varepsilon}(\boldsymbol u),
\qquad
\tilde{\boldsymbol E}
=
\dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u}).
$$

代码中使用的第二不变量和正则化 Glen 黏度为

$$
\varepsilon_{II}
=
\frac12\boldsymbol E:\boldsymbol E
=
\frac12
\left(
E_{xx}^2+E_{zz}^2+2E_{xz}^2
\right),
$$

$$
\eta
=
\frac12 A^{-1/n}
s^\alpha,
\qquad
s
=
\varepsilon_{II}+\varepsilon_{\rm reg}^2,
\qquad
\alpha
=
\frac{1-n}{2n}.
$$

黏性应力为

$$
\boldsymbol\sigma_{\rm visc}
=
2\eta\boldsymbol E.
$$

对 $\tilde{\boldsymbol E}$ 做一阶线性化，得到

$$
\mathbb C(\boldsymbol u):\tilde{\boldsymbol E}
=
2\eta\tilde{\boldsymbol E}
+
2\eta
\frac{\alpha}{s}
\left(
\boldsymbol E:\tilde{\boldsymbol E}
\right)
\boldsymbol E.
$$

因此四阶张量分量可以写为

$$
C_{ijkl}(\boldsymbol u)
=
2\eta\,\delta_{ik}\delta_{jl}
+
2\eta
\frac{\alpha}{s}
E_{ij}E_{kl}.
$$

这里第一项是冻结黏度项；第二项来自黏度 $\eta$ 对当前应变率的导数。若只保留第一项，就对应 Picard 线性化；完整两项一起才是一致切线。

在代码中并没有显式保存四阶张量 $C_{ijkl}$。`NonlinearStokesP2P1.m` 中的 `assembleviscoustangent(uk)` 直接装配矩阵条目。对测试函数方向 $\phi_a$ 和试探函数方向 $\phi_b$，代码先计算

$$
\dot{\boldsymbol\varepsilon}_a
=
\dot{\boldsymbol\varepsilon}(\phi_a),
\qquad
\dot{\boldsymbol\varepsilon}_b
=
\dot{\boldsymbol\varepsilon}(\phi_b).
$$

代码中的

$$
\texttt{strainDot}
=
\dot{\boldsymbol\varepsilon}_a:
\dot{\boldsymbol\varepsilon}_b,
$$

$$
\texttt{stateDotA}
=
\boldsymbol E:
\dot{\boldsymbol\varepsilon}_a,
\qquad
\texttt{stateDotB}
=
\boldsymbol E:
\dot{\boldsymbol\varepsilon}_b.
$$

于是黏性一致切线矩阵的局部条目为

$$
K_{ab}^{\rm tan}
=
\int_\Omega
2\eta
\left[
\dot{\boldsymbol\varepsilon}_a:
\dot{\boldsymbol\varepsilon}_b
+
\frac{\alpha}{s}
\left(
\boldsymbol E:\dot{\boldsymbol\varepsilon}_a
\right)
\left(
\boldsymbol E:\dot{\boldsymbol\varepsilon}_b
\right)
\right]\,dx.
$$

这正是

$$
K_{ab}^{\rm tan}
=
\int_\Omega
\dot{\boldsymbol\varepsilon}_a:
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_b
\,dx.
$$

对应的 MATLAB 代码是

```matlab
stateDotA = exx.*aExx+ezz.*aEzz+2*exz.*aExz;
stateDotB = exx.*bExx+ezz.*bEzz+2*exz.*bExz;
strainDot = aExx.*bExx+aEzz.*bEzz+2*aExz.*bExz;

kab = 2*w(q)*area.*eta.*...
    (strainDot+exponent./strainRegularized.*...
     stateDotA.*stateDotB);
```

其中 `exponent` 就是 $\alpha$，`strainRegularized` 就是 $s$。最后这些局部条目装配成黏性一致切线矩阵 `Kt`，再与底部滑移切线矩阵 `Kbt`、散度矩阵和约束矩阵组合为程序中的

$$
\texttt{eqn.tangent}.
$$

### 底部滑移线性化

因为本文取 $m=1$，底部摩擦项是

$$
\beta\boldsymbol u_t.
$$

对状态方向 $\tilde{\boldsymbol u}$，有

$$
D_{\boldsymbol u}
\left(
\beta\boldsymbol u_t
\right)
[\tilde{\boldsymbol u}]
=
\beta\tilde{\boldsymbol u}_t.
$$

对参数方向 $\delta\beta$，有

$$
D_\beta
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
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p]
=
-
\mathcal R_\beta[\delta\beta],
$$

得到强形式

$$
\left\{
\begin{aligned}
\nabla\cdot\tilde{\boldsymbol u} &=0
&&\text{in }\Omega,\\
-\nabla\cdot
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
&=\boldsymbol 0
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol u}\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
+
\beta\tilde{\boldsymbol u}_t
&=
-\delta\beta\boldsymbol u_t
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol u},\
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

右端只出现在底部，是因为 $\beta$ 只出现在底部滑移边界条件中。

## 伴随问题

如果直接计算梯度 $\delta J/\delta\beta$，需要知道每个参数方向 $\delta\beta$ 对应的状态扰动 $\tilde{\boldsymbol u}$。伴随方法的目的，是用一个伴随问题消去 $\tilde{\boldsymbol u}$，从而直接得到关于 $\delta\beta$ 的梯度表达。

令增量正问题的线性算子为

$$
\mathcal A(\tilde{\boldsymbol u},\tilde p)
=
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p].
$$

伴随算子 $\mathcal A^*$ 由弱形式内积定义：

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

这就是“伴随方程使用线性化算子的转置”的含义。有限维类比是：增量正问题使用 $A$，伴随问题使用 $A^T$。

为了抵消目标函数变分

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)\tilde u_x\,ds,
$$

伴随顶部边界条件取为

$$
\delta\boldsymbol\sigma^*[\boldsymbol v,r]\boldsymbol n
=
-
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\qquad\text{on }\Gamma_t.
$$

其中负号来自把目标函数变分移到伴随方程边界项的另一侧。若采用不同的拉格朗日函数符号约定，伴随变量可能整体差一个负号，但梯度公式会相应保持一致。

### 伴随线性化应力

若

$$
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u})
-
\tilde p\boldsymbol I,
$$

则伴随线性化应力定义为

$$
\delta\boldsymbol\sigma^*[\boldsymbol v,r]
=
\mathbb C(\boldsymbol u)^T:
\dot{\boldsymbol\varepsilon}(\boldsymbol v)
-
r\boldsymbol I.
$$

这里 $\mathbb C(\boldsymbol u)^T$ 是相对于张量内积的转置，即

$$
\left(
\mathbb C:\boldsymbol E
\right):\boldsymbol F
=
\boldsymbol E:
\left(
\mathbb C^T:\boldsymbol F
\right).
$$

对于由标量黏性势导出的 Glen 型黏性切线，$\mathbb C$ 通常是对称的，所以 $\mathbb C^T=\mathbb C$。这时体内的伴随线性化应力在形式上和增量正问题中的线性化应力相同，只是方向变量换成 $(\boldsymbol v,r)$。

### 伴随方程

伴随问题为求 $(\boldsymbol v,r)$，满足

$$
\left\{
\begin{aligned}
\nabla\cdot\boldsymbol v &=0
&&\text{in }\Omega,\\
-\nabla\cdot
\delta\boldsymbol\sigma^*[\boldsymbol v,r]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma^*[\boldsymbol v,r]\boldsymbol n
&=
-\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
&&\text{on }\Gamma_t,\\
\boldsymbol v\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma^*[\boldsymbol v,r]\boldsymbol n
+
\beta\boldsymbol v_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol v,\
\delta\boldsymbol\sigma^*[\boldsymbol v,r]\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

底部项中出现 $\beta\boldsymbol v_t$，是因为 $m=1$ 时底部滑移线性化算子为 $\beta\boldsymbol T$，而 $\boldsymbol T$ 是自伴随投影。

## 梯度公式

由增量正问题和伴随问题相消，可以把目标函数对参数方向的变化写成底部边界积分：

$$
\delta J[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

如果使用 $q=\log\beta$，则 $\delta\beta=\beta\delta q$，于是

$$
\delta J[\delta q]
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
\nabla\cdot\tilde{\boldsymbol v} &=0
&&\text{in }\Omega,\\
-\nabla\cdot
\delta\boldsymbol\sigma^*[\tilde{\boldsymbol v},\tilde r]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma^*[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
&=
-\tilde u_x\boldsymbol e_x
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol v}\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma^*[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
+
\beta\tilde{\boldsymbol v}_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol v},\
\delta\boldsymbol\sigma^*[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
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

## 四个方程的逻辑关系

四个方程可以压缩成下面的抽象关系。

正问题：

$$
\mathcal R(\boldsymbol u,p,\beta)=0.
$$

增量正问题：

$$
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p]
=
-
\mathcal R_\beta[\delta\beta].
$$

伴随问题：

$$
\mathcal R_{(\boldsymbol u,p)}^*
[\boldsymbol v,r]
=
-
J_{\boldsymbol u}.
$$

增量伴随问题：

$$
\mathcal R_{(\boldsymbol u,p)}^*
[\tilde{\boldsymbol v},\tilde r]
=
-
J_{\boldsymbol u\boldsymbol u}
[\tilde{\boldsymbol u}].
$$

这里 $J_{\boldsymbol u}$ 是目标函数对状态的导数，$J_{\boldsymbol u\boldsymbol u}[\tilde{\boldsymbol u}]$ 是目标函数二阶导数作用在增量状态方向上。

对于顶部水平速度观测，有

$$
J_{\boldsymbol u}
=
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\quad\text{on }\Gamma_t,
$$

以及

$$
J_{\boldsymbol u\boldsymbol u}
[\tilde{\boldsymbol u}]
=
\tilde u_x\boldsymbol e_x
\quad\text{on }\Gamma_t.
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
\mathcal R_\beta[\delta\beta].
$$

然后引入伴随变量 $(\boldsymbol v,r)$，令它满足线性化算子的转置问题：

$$
\mathcal A^*(\boldsymbol v,r)
=
-
J_{\boldsymbol u}.
$$

这样目标函数变分中的状态扰动 $\tilde{\boldsymbol u}$ 可以被消去，最终梯度只表现为底部边界上的参数项：

$$
\delta J[\delta\beta]
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

| 符号                                                            | 含义                 |
| :---:                                                           | :---:                |
| $\Omega$                                                        | 冰体区域             |
| $\Gamma_t$                                                      | 上表面               |
| $\Gamma_b$                                                      | 底部                 |
| $\Gamma_p$                                                      | 周期边界             |
| $\boldsymbol n$                                                 | 外法向               |
| $\boldsymbol T=\boldsymbol I-\boldsymbol n\otimes\boldsymbol n$ | 切向投影             |
| $\boldsymbol u$                                                 | 正问题速度           |
| $p$                                                             | 正问题压力           |
| $\boldsymbol v$                                                 | 伴随速度             |
| $r$                                                             | 伴随压力             |
| $\tilde{\boldsymbol u}$                                         | 增量正问题速度方向   |
| $\tilde p$                                                      | 增量正问题压力方向   |
| $\tilde{\boldsymbol v}$                                         | 增量伴随速度         |
| $\tilde r$                                                      | 增量伴随压力         |
| $\beta$                                                         | 底部摩擦系数         |
| $q=\log\beta$                                                   | 对数参数             |
| $\delta\beta$                                                   | 摩擦系数扰动         |
| $\delta q$                                                      | 对数参数扰动         |
| $\dot{\boldsymbol\varepsilon}(\boldsymbol u)$                   | 应变率张量           |
| $\eta(\boldsymbol u)$                                           | 有效黏度             |
| $\boldsymbol\sigma$                                             | 正问题应力           |
| $\delta\boldsymbol\sigma$                                       | 应力的一阶线性化     |
| $\delta\boldsymbol\sigma^*$                                     | 线性化应力算子的伴随 |
| $u_{\rm obs}$                                                   | 顶部水平速度观测     |
| $J$                                                             | 反问题目标函数       |
| $\mathbb C(\boldsymbol u)$                                      | 黏性一致切线张量     |
