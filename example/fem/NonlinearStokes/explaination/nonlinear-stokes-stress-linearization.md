# 非线性 Stokes 反问题中的四个方程

本文从数学角度整理非线性 Stokes 冰流反问题中常见的四类方程：

1. 正问题；
2. 伴随问题；
3. 增量正问题；
4. 增量伴随问题。

核心思想是：正问题给定底部摩擦系数 $\beta$，求速度压力 $(\boldsymbol u,p)$；反问题根据顶部观测速度恢复 $\beta$。伴随问题来自目标函数对正问题约束的拉格朗日乘子；增量正问题来自正问题对参数扰动的一阶线性化；增量伴随问题来自伴随问题或梯度对参数方向的进一步线性化，通常用于构造 Gauss--Newton Hessian 作用。

全文只讨论连续数学形式，不涉及离散矩阵或代码实现。

## 几何区域和边界

设冰体区域为 $\Omega\subset\mathbb R^2$，边界分成三部分：

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

对任意速度 $\boldsymbol u$，它的底部切向分量为

$$
\boldsymbol u_t
=
\boldsymbol T\boldsymbol u.
$$

底部不可穿透条件是

$$
\boldsymbol u\cdot\boldsymbol n=0
\qquad\text{on }\Gamma_b.
$$

底部滑移律只作用在切向方向。

## 应变率、黏度和应力

速度应变率张量为

$$
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
=
\frac12
\left(
\nabla\boldsymbol u+\nabla\boldsymbol u^T
\right).
$$

为了简化记号，记

$$
\boldsymbol E
=
\dot{\boldsymbol\varepsilon}(\boldsymbol u).
$$

非线性 Stokes 模型中的应力为

$$
\boldsymbol\sigma
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
-
p\boldsymbol I.
$$

这里：

- $\eta(\boldsymbol u)$ 是有效黏度；
- $\eta$ 依赖于应变率不变量，因此间接依赖 $\boldsymbol u$；
- $p$ 是压力；
- $-\nabla\cdot\boldsymbol\sigma$ 是动量方程中的内力项。

例如 Glen 型流律可以抽象写成

$$
\eta(\boldsymbol u)
=
\eta\left(
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
\right).
$$

本文推导不需要固定 $\eta$ 的具体表达式，只需要知道它依赖 $\boldsymbol u$。

## 底部滑移律

底部摩擦系数为 $\beta$。设滑移指数为 $m$。底部切向摩擦力写为

$$
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t.
$$

因此底部滑移边界条件是

$$
\boldsymbol T\boldsymbol\sigma\boldsymbol n
+
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
=
\boldsymbol 0
\qquad\text{on }\Gamma_b.
$$

其中第一项 $\boldsymbol T\boldsymbol\sigma\boldsymbol n$ 是底部牵引的切向分量，第二项是底部摩擦力。

## 正问题

给定底部摩擦系数 $\beta$，正问题是求速度压力 $(\boldsymbol u,p)$，满足

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
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol u,\ \boldsymbol\sigma\boldsymbol n
&\text{ satisfy periodicity}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

这里

$$
\boldsymbol\sigma
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}(\boldsymbol u)
-
p\boldsymbol I.
$$

正问题可以看成一个非线性约束

$$
\mathcal R(\boldsymbol u,p,\beta)=0.
$$

反问题中的所有伴随和增量方程，都是从这个约束的一阶线性化和转置得到的。

## 目标函数

假设顶部有水平速度观测 $u_{\rm obs}$。目标函数取为

$$
J(\boldsymbol u)
=
\frac12
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)^2\,ds.
$$

其中

$$
u_x=\boldsymbol u\cdot\boldsymbol e_x.
$$

目标函数的一阶变分是

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)
\tilde u_x\,ds.
$$

因为

$$
\tilde u_x
=
\tilde{\boldsymbol u}\cdot\boldsymbol e_x,
$$

所以也可以写成

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

这个顶部边界项就是伴随问题中顶部边界条件的来源。

## 应力的一阶线性化

令状态扰动方向为

$$
\tilde{\boldsymbol u},
\qquad
\tilde p.
$$

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

由于

$$
\boldsymbol\sigma(\boldsymbol u_\epsilon,p_\epsilon)
=
2\eta(\boldsymbol u_\epsilon)
\dot{\boldsymbol\varepsilon}(\boldsymbol u_\epsilon)
-
p_\epsilon\boldsymbol I,
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

这三个项分别来自：

- 应变率变化；
- 黏度随速度变化；
- 压力变化。

如果把黏性部分写成四阶切线张量形式，则有

$$
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}(\tilde{\boldsymbol u})
-
\tilde p\boldsymbol I.
$$

其中 $\mathbb C(\boldsymbol u)$ 是非线性黏性应力对 $\dot{\boldsymbol\varepsilon}(\boldsymbol u)$ 的导数。这个张量就是连续层面的一致切线。

## 底部滑移律的一阶线性化

底部摩擦项为

$$
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t.
$$

对速度扰动 $\tilde{\boldsymbol u}$，因为

$$
\tilde{\boldsymbol u}_t
=
\boldsymbol T\tilde{\boldsymbol u},
$$

有

$$
\begin{aligned}
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)
[\tilde{\boldsymbol u}]
&=
\beta|\boldsymbol u_t|^{m-1}\tilde{\boldsymbol u}_t
\\
&\quad
+
\beta(m-1)|\boldsymbol u_t|^{m-3}
\left(\boldsymbol u_t\cdot\tilde{\boldsymbol u}_t\right)
\boldsymbol u_t.
\end{aligned}
$$

对参数扰动 $\delta\beta$，有

$$
D_\beta
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)
[\delta\beta]
=
\delta\beta
|\boldsymbol u_t|^{m-1}\boldsymbol u_t.
$$

如果使用对数参数

$$
q=\log\beta,
$$

那么

$$
\delta\beta
=
\beta\,\delta q.
$$

## 增量正问题

增量正问题描述：当参数 $\beta$ 沿方向 $\delta\beta$ 改变时，正问题解 $(\boldsymbol u,p)$ 的一阶变化 $(\tilde{\boldsymbol u},\tilde p)$ 满足什么方程。

从正问题残差

$$
\mathcal R(\boldsymbol u,p,\beta)=0
$$

出发，对 $(\boldsymbol u,p,\beta)$ 做一阶变分：

$$
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p]
+
\mathcal R_\beta[\delta\beta]
=
0.
$$

因此

$$
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p]
=
-
\mathcal R_\beta[\delta\beta].
$$

写成强形式就是

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
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)
[\tilde{\boldsymbol u}]
&=
-\delta\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol u},\
\delta\boldsymbol\sigma[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
&\text{ satisfy periodicity}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

右端只出现在底部，是因为参数 $\beta$ 只出现在底部滑移边界条件中。

## 伴随问题的来源

为了避免直接计算 $\delta J/\delta\beta$ 中的状态导数，引入伴随变量 $(\boldsymbol v,r)$。

令增量正问题的线性算子为

$$
\mathcal A(\tilde{\boldsymbol u},\tilde p)
=
\mathcal R_{(\boldsymbol u,p)}
[\tilde{\boldsymbol u},\tilde p].
$$

伴随算子 $\mathcal A^*$ 由如下恒等式定义：

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

伴随变量的作用是让目标函数变分

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
$$

被线性化正问题的伴随边界项抵消。由于

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
\left(u_x-u_{\rm obs}\right)\tilde u_x\,ds,
$$

伴随问题在顶部产生边界条件

$$
\delta\boldsymbol\sigma^*[\boldsymbol v,r]\boldsymbol n
=
-
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\qquad\text{on }\Gamma_t.
$$

这里负号来自把目标函数变分移到伴随方程边界项的另一侧。若采用不同的拉格朗日函数符号约定，伴随变量整体可能差一个负号，但最终梯度表达一致。

## 伴随线性化应力

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

其中 $\mathbb C(\boldsymbol u)^T$ 表示相对于张量内积的转置，即

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

对于由标量黏性势导出的 Glen 型黏性切线，$\mathbb C$ 通常是对称的，因此

$$
\mathbb C^T=\mathbb C.
$$

这时体内的伴随线性化应力在形式上和增量正问题中的线性化应力相同，只是把方向 $(\tilde{\boldsymbol u},\tilde p)$ 换成 $(\boldsymbol v,r)$。

## 伴随问题

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
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)^*
[\boldsymbol v]
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol v,\
\delta\boldsymbol\sigma^*[\boldsymbol v,r]\boldsymbol n
&\text{ satisfy periodicity}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

其中 $D_{\boldsymbol u}(\cdot)^*$ 是底部滑移线性化算子的伴随。

对滑移律

$$
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t,
$$

其速度导数对应的切向矩阵为

$$
\beta|\boldsymbol u_t|^{m-1}\boldsymbol T
+
\beta(m-1)|\boldsymbol u_t|^{m-3}
\boldsymbol u_t\otimes\boldsymbol u_t.
$$

这个矩阵关于欧氏内积是对称的，所以在这种滑移律下

$$
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)^*
[\boldsymbol v]
=
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)
[\boldsymbol v].
$$

因此伴随底部边界条件与增量正问题左端在形式上相同。

## 梯度公式

目标函数对参数方向 $\delta\beta$ 的一阶变化可以通过伴随变量写成底部边界积分。

由增量正问题和伴随问题相消，得到

$$
\delta J[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta
|\boldsymbol u_t|^{m-1}
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds,
$$

其中

$$
\boldsymbol v_t
=
\boldsymbol T\boldsymbol v.
$$

如果使用 $q=\log\beta$，则 $\delta\beta=\beta\delta q$，于是

$$
\delta J[\delta q]
=
\int_{\Gamma_b}
\beta\delta q
|\boldsymbol u_t|^{m-1}
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

所以相对于 $q$ 的连续梯度密度为

$$
g_q
=
\beta
|\boldsymbol u_t|^{m-1}
\boldsymbol u_t\cdot\boldsymbol v_t.
$$

符号可能随伴随变量的符号约定整体改变；只要伴随方程和梯度公式使用同一约定，下降方向是一致的。

## 增量伴随问题的来源

增量伴随问题用于计算梯度沿参数方向 $\delta\beta$ 的变化，也就是 Hessian-vector product。

给定参数方向 $\delta\beta$，先由增量正问题得到状态方向

$$
(\tilde{\boldsymbol u},\tilde p).
$$

目标函数二阶变分中，观测项对状态方向的导数是

$$
\delta
\left[
\left(u_x-u_{\rm obs}\right)\boldsymbol e_x
\right]
=
\tilde u_x\boldsymbol e_x.
$$

因此增量伴随问题顶部边界条件由 $\tilde u_x$ 驱动。

Gauss--Newton 近似忽略正问题二阶导数和伴随方程中由非线性算子二阶导数带来的项，只保留观测算子的二阶项。于是增量伴随问题为求 $(\tilde{\boldsymbol v},\tilde r)$，满足

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
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)^*
[\tilde{\boldsymbol v}]
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol v},\
\delta\boldsymbol\sigma^*[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
&\text{ satisfy periodicity}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

它和伴随问题的区别只在顶部边界源项：

$$
-(u_x-u_{\rm obs})\boldsymbol e_x
\quad\longrightarrow\quad
-\tilde u_x\boldsymbol e_x.
$$

前者由当前观测误差驱动，后者由增量状态在观测量上的变化驱动。

## 四个方程的逻辑关系

四个方程可以按如下逻辑理解：

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

## 符号表

| 符号 | 含义 |
|---|---|
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
| $\dot{\boldsymbol\varepsilon}(\boldsymbol u)$ | 应变率张量 |
| $\eta(\boldsymbol u)$ | 有效黏度 |
| $\boldsymbol\sigma$ | 正问题应力 |
| $\delta\boldsymbol\sigma$ | 应力的一阶线性化 |
| $\delta\boldsymbol\sigma^*$ | 线性化应力算子的伴随 |
| $u_{\rm obs}$ | 顶部水平速度观测 |
| $J$ | 反问题目标函数 |
| $\mathbb C(\boldsymbol u)$ | 黏性一致切线张量 |
