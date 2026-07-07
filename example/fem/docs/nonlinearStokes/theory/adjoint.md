# 伴随问题的推导

本文补充说明 [equations.md](equations.md) 中伴随问题的来源。目标是把增量正问题中的各个线性化算子逐项取伴随，并解释连续层面的伴随方程、边界条件和梯度公式。

## 1. 从增量正问题开始

在当前正问题解 $(\boldsymbol u,p)$ 和参数 $\beta$ 固定时，状态方向 $(\tilde{\boldsymbol u},\tilde p)$ 的线性化应力为

$$
\delta\boldsymbol\sigma(\boldsymbol u,p)
[\tilde{\boldsymbol u},\tilde p]
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
-
\tilde p\boldsymbol I.
$$

下面简记

$$
\boldsymbol\tau
=
\boldsymbol\tau(\tilde{\boldsymbol u},\tilde p)
=
\delta\boldsymbol\sigma(\boldsymbol u,p)
[\tilde{\boldsymbol u},\tilde p].
$$

增量正问题的状态算子可以按部件理解为

$$
\mathcal A(\tilde{\boldsymbol u},\tilde p)
=
\left(
\begin{array}{c}
-\nabla\cdot\boldsymbol\tau\\
\nabla\cdot\tilde{\boldsymbol u}\\
\boldsymbol\tau\boldsymbol n|_{\Gamma_t}\\
\tilde{\boldsymbol u}\cdot\boldsymbol n|_{\Gamma_b}\\
\left(
\boldsymbol T\boldsymbol\tau\boldsymbol n
+\beta\tilde{\boldsymbol u}_t
\right)|_{\Gamma_b}\\
\text{periodic jump}
\end{array}
\right).
$$

参数方向 $\delta\beta$ 只出现在底部滑移条件中：

$$
\boldsymbol T\boldsymbol\tau\boldsymbol n
+\beta\tilde{\boldsymbol u}_t
=
-\delta\beta\,\boldsymbol u_t
\qquad\text{on }\Gamma_b.
$$

目标函数的一阶变分是顶部项

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
(u_x-u_{\rm obs})\boldsymbol e_x
\cdot\tilde{\boldsymbol u}\,ds.
$$

伴随推导要做的事就是：用一个伴随变量 $(\boldsymbol v,r)$ 把这个顶部的 $\tilde{\boldsymbol u}$ 消掉，最后只留下底部的 $\delta\beta$ 项。

## 2. 固定弱形式配对并引入 $\mathcal A^*$

伴随算子不是只由强形式微分表达式决定，还依赖采用什么内积和符号约定。本文把增量正问题的完整状态残差和伴随变量放在同一个配对中：

$$
\begin{aligned}
I
&=
\left\langle
\mathcal A(\tilde{\boldsymbol u},\tilde p),
(\boldsymbol v,r)
\right\rangle\\
&=
\int_\Omega
\left(-\nabla\cdot\boldsymbol\tau\right)
\cdot\boldsymbol v\,dx
-
\int_\Omega
r\,\nabla\cdot\tilde{\boldsymbol u}\,dx\\
&=: I_1+I_2.
\end{aligned}
$$

这里 $I$ 不是新的物理量，而是定义伴随算子的中间表达式。分部积分以后，要把它整理成

$$
I
=
\left\langle
(\tilde{\boldsymbol u},\tilde p),
\mathcal A^*(\boldsymbol v,r)
\right\rangle
+
\text{boundary terms}.
$$

因此，$\mathcal A^*(\boldsymbol v,r)$ 就是分部积分后乘在任意方向 $(\tilde{\boldsymbol u},\tilde p)$ 前面的体内系数。后面会得到

$$
\mathcal A^*(\boldsymbol v,r)
=
\left(
\begin{array}{c}
-\nabla\cdot\boldsymbol\tau^*(\boldsymbol v,r)\\
-\nabla\cdot\boldsymbol v
\end{array}
\right),
$$

其中

$$
\boldsymbol\tau^*(\boldsymbol v,r)
=
\delta\boldsymbol\sigma^*(\boldsymbol u,p)[\boldsymbol v,r]
=
\mathbb C(\boldsymbol u)^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
-
r\boldsymbol I.
$$

之所以要把动量方程和质量方程放在同一个 $I$ 里，是因为不可压 Stokes 是速度和压力耦合的 saddle-point 系统。增量正问题的状态未知量是一个整体 $(\tilde{\boldsymbol u},\tilde p)$，对应的线性化残差也必须作为一个整体。伴随变量同样是一个整体 $(\boldsymbol v,r)$。

第一项是增量动量方程和伴随速度 $\boldsymbol v$ 的配对。第二项是增量不可压缩方程和伴随压力 $r$ 的配对，前面的负号是本文的符号约定。如果把第二项改成 $+\int_\Omega r\nabla\cdot\tilde{\boldsymbol u}\,dx$，伴随压力的符号会相应改变。只要伴随方程和梯度公式使用同一套约定，结果是一致的。

## 3. 计算配对并读出体内伴随算子

先处理动量配对 $I_1$。分部积分得

$$
\int_\Omega
\left(-\nabla\cdot\boldsymbol\tau\right)
\cdot\boldsymbol v\,dx
=
\int_\Omega
\boldsymbol\tau:\nabla\boldsymbol v\,dx
-
\int_{\partial\Omega}
(\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v\,ds.
$$

由于 $\boldsymbol\tau$ 是对称应力，

$$
\boldsymbol\tau:\nabla\boldsymbol v
=
\boldsymbol\tau:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}.
$$

代入 $\boldsymbol\tau=\mathbb C:\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}-\tilde p\boldsymbol I$，并使用四阶张量转置定义

$$
(\mathbb C:\boldsymbol A):\boldsymbol B
=
\boldsymbol A:(\mathbb C^T:\boldsymbol B),
$$

得到

$$
I_1
=
\int_\Omega
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
:
\left(
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
\right)\,dx
-
\int_\Omega
\tilde p\,\nabla\cdot\boldsymbol v\,dx
-
\int_{\partial\Omega}
(\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v\,ds.
$$

现在速度方向还以 $\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$ 的形式出现。为了读出强形式伴随动量方程，需要把这个导数从 $\tilde{\boldsymbol u}$ 上移走。令

$$
\boldsymbol M
=
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}.
$$

在本文的 Stokes 应力线性化中，$\boldsymbol M$ 是对称张量，所以

$$
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}:\boldsymbol M
=
\nabla\tilde{\boldsymbol u}:\boldsymbol M.
$$

由张量分部积分公式，

$$
\int_\Omega
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
:
\left(
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
\right)\,dx
=
-\int_\Omega
\nabla\cdot
\left(
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
\right)
\cdot\tilde{\boldsymbol u}\,dx
+
\int_{\partial\Omega}
\left(
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
\right)\boldsymbol n
\cdot\tilde{\boldsymbol u}\,ds.
$$

于是

$$
\begin{aligned}
I_1
&=
-\int_\Omega
\nabla\cdot
\left(
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
\right)
\cdot\tilde{\boldsymbol u}\,dx
-
\int_\Omega
\tilde p\,\nabla\cdot\boldsymbol v\,dx\\
&\quad+
\int_{\partial\Omega}
\left(
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
\right)\boldsymbol n
\cdot\tilde{\boldsymbol u}\,ds
-
\int_{\partial\Omega}
(\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v\,ds.
\end{aligned}
$$

再处理不可压缩配对 $I_2$：

$$
I_2
=
-\int_\Omega
r\,\nabla\cdot\tilde{\boldsymbol u}\,dx
=
\int_\Omega
\nabla r\cdot\tilde{\boldsymbol u}\,dx
-
\int_{\partial\Omega}
r\,\tilde{\boldsymbol u}\cdot\boldsymbol n\,ds.
$$

定义

$$
\boldsymbol\tau^*(\boldsymbol v,r)
=
\mathbb C^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
-
r\boldsymbol I.
$$

合并 $I=I_1+I_2$，得到

$$
\begin{aligned}
I
&=
\int_\Omega
\left[-\nabla\cdot\boldsymbol\tau^*(\boldsymbol v,r)\right]
\cdot\tilde{\boldsymbol u}\,dx
+
\int_\Omega
\tilde p\left[-\nabla\cdot\boldsymbol v\right]\,dx\\
&\quad+
\int_{\partial\Omega}
\left[
\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
\cdot\tilde{\boldsymbol u}
-
\boldsymbol\tau\boldsymbol n
\cdot\boldsymbol v
\right]\,ds.
\end{aligned}
$$

这正是第 2 节中

$$
I
=
\left\langle
(\tilde{\boldsymbol u},\tilde p),
\mathcal A^*(\boldsymbol v,r)
\right\rangle
+
\text{boundary terms}
$$

的具体形式。因此体内伴随算子为

$$
\mathcal A^*(\boldsymbol v,r)
=
\left(
\begin{array}{c}
-\nabla\cdot\boldsymbol\tau^*(\boldsymbol v,r)\\
-\nabla\cdot\boldsymbol v
\end{array}
\right).
$$

伴随法要把目标函数变分中的状态增量
$\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]$ 消去。分部积分后，体内所有状态方向项已经集中为

$$
\left\langle
(\tilde{\boldsymbol u},\tilde p),
\mathcal A^*(\boldsymbol v,r)
\right\rangle
$$

要让最终公式不含未知的 $(\tilde{\boldsymbol u},\tilde p)$，需要这项对任意状态方向都为零。因此取

$$
\mathcal A^*(\boldsymbol v,r)=0
\qquad\text{in }\Omega.
$$

代入上面的 $\mathcal A^*$，体内伴随方程为

$$
-\nabla\cdot\boldsymbol\tau^*(\boldsymbol v,r)
=
\boldsymbol 0,
\qquad
\nabla\cdot\boldsymbol v=0
\qquad\text{in }\Omega.
$$

对本文的 Glen 型黏性切线，$\mathbb C^T=\mathbb C$，所以伴随应力和增量正问题的线性化应力形式相同，只是方向变量换成 $(\boldsymbol v,r)$。

## 4. 边界项的统一形式

经过两次分部积分后，除了体积分，还剩边界项

$$
-\int_{\partial\Omega}
(\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v\,ds
+
\int_{\partial\Omega}
(\boldsymbol\tau^*\boldsymbol n)\cdot\tilde{\boldsymbol u}\,ds.
$$

伴随边界条件就是从这个式子在不同边界上的消去要求读出来的。

## 5. 顶部边界：观测误差成为伴随牵引

在顶部，增量正问题满足

$$
\boldsymbol\tau\boldsymbol n=\boldsymbol 0
\qquad\text{on }\Gamma_t.
$$

因此顶部边界项只剩

$$
\int_{\Gamma_t}
(\boldsymbol\tau^*\boldsymbol n)
\cdot\tilde{\boldsymbol u}\,ds.
$$

约化目标函数的状态变分为

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_t}
(u_x-u_{\rm obs})\boldsymbol e_x
\cdot\tilde{\boldsymbol u}\,ds
$$

伴随变量的选择要使这个顶部状态增量项被边界项抵消，即

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
+
\int_{\Gamma_t}
(\boldsymbol\tau^*\boldsymbol n)
\cdot\tilde{\boldsymbol u}\,ds
=0
$$

对任意 $\tilde{\boldsymbol u}$ 成立。于是得到顶部伴随牵引

$$
\boxed{
\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
=
-(u_x-u_{\rm obs})\boldsymbol e_x
\qquad\text{on }\Gamma_t.
}
$$

这里的负号来自“顶部边界项抵消目标函数变分”的符号约定。

## 6. 底部边界：滑移算子的伴随和梯度

底部增量速度和伴随速度取在同一个不可穿透空间中：

$$
\tilde{\boldsymbol u}\cdot\boldsymbol n=0
\qquad
\boldsymbol v\cdot\boldsymbol n=0
\qquad\text{on }\Gamma_b.
$$

因此底部只需要比较切向分量。由于 $\boldsymbol v=\boldsymbol v_t$，有

$$
(\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v
=
(\boldsymbol T\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v_t.
$$

增量滑移条件给出

$$
\boldsymbol T\boldsymbol\tau\boldsymbol n
=
-\beta\tilde{\boldsymbol u}_t
-\delta\beta\,\boldsymbol u_t.
$$

代入底部边界项

$$
-\int_{\Gamma_b}
(\boldsymbol\tau\boldsymbol n)\cdot\boldsymbol v\,ds
+
\int_{\Gamma_b}
(\boldsymbol\tau^*\boldsymbol n)\cdot\tilde{\boldsymbol u}\,ds,
$$

得到

$$
\int_{\Gamma_b}
\left(
\boldsymbol T\boldsymbol\tau^*\boldsymbol n
+\beta\boldsymbol v_t
\right)
\cdot\tilde{\boldsymbol u}_t\,ds
+
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t\,ds.
$$

为了消去任意的底部切向状态方向 $\tilde{\boldsymbol u}_t$，伴随底部滑移条件必须是

$$
\boxed{
\boldsymbol T\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
+\beta\boldsymbol v_t
=
\boldsymbol 0
\qquad\text{on }\Gamma_b.
}
$$

消去 $\tilde{\boldsymbol u}_t$ 后，底部只剩参数项

$$
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t\,ds.
$$

下面说明它为什么就是方向导数。对参数求导时，目标函数应理解为约化目标函数

$$
J(\beta):=J(\boldsymbol u(\beta)).
$$

链式法则给出

$$
\delta J(\beta)[\delta\beta]
=
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}].
$$

其中 $(\tilde{\boldsymbol u},\tilde p)$ 是 $\delta\beta$ 驱动的增量正问题解。增量正问题可写成

$$
\mathcal A(\tilde{\boldsymbol u},\tilde p)
+
\mathcal R_\beta[\delta\beta]
=0,
$$

其中底部滑移条件中的参数残差是

$$
\mathcal R_\beta[\delta\beta]
=
\delta\beta\,\boldsymbol u_t
\qquad\text{on }\Gamma_b.
$$

把增量方程与伴随变量配对并分部积分后，伴随体方程消去体积分，顶部伴随牵引给出 $-\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]$，底部伴随滑移条件消去 $\tilde{\boldsymbol u}_t$。因此配对恒等式化为

$$
0
=
-\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
+
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t\,ds.
$$

于是

$$
\delta J(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t\,ds.
$$

再由链式法则得到梯度公式

$$
\delta J(\beta)[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta\,
\boldsymbol u_t\cdot\boldsymbol v_t\,ds.
$$

若使用 $q=\log\beta$，则 $\delta\beta=\beta\delta q$，所以

$$
\delta J(q)[\delta q]
=
\int_{\Gamma_b}
\beta\,\delta q\,
\boldsymbol u_t\cdot\boldsymbol v_t\,ds,
\qquad
g_q
=
\beta\boldsymbol u_t\cdot\boldsymbol v_t.
$$

## 7. 周期边界

在周期边界 $\Gamma_p$ 上，增量正问题要求速度和牵引周期匹配。分部积分产生的两侧边界项方向相反。

若伴随变量满足

$$
\boldsymbol v,\quad
\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
\quad\text{periodic on }\Gamma_p,
$$

则周期边界上的两侧积分相互抵消。因此周期条件的伴随仍是周期条件。

## 8. 伴随强形式汇总

伴随问题为求 $(\boldsymbol v,r)$，使

$$
\left\{
\begin{aligned}
\nabla\cdot\boldsymbol v &=0
&&\text{in }\Omega,\\
-\nabla\cdot
\boldsymbol\tau^*(\boldsymbol v,r)
&=\boldsymbol 0
&&\text{in }\Omega,\\
\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
&=
-(u_x-u_{\rm obs})\boldsymbol e_x
&&\text{on }\Gamma_t,\\
\boldsymbol v\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
+\beta\boldsymbol v_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\boldsymbol v,\quad
\boldsymbol\tau^*(\boldsymbol v,r)\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

其中

$$
\boldsymbol\tau^*(\boldsymbol v,r)
=
\mathbb C(\boldsymbol u)^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
-
r\boldsymbol I.
$$

若 $\mathbb C^T=\mathbb C$，则

$$
\boldsymbol\tau^*(\boldsymbol v,r)
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
-
r\boldsymbol I.
$$

## 9. 算子伴随对照

| 增量正问题中的项 | 伴随计算 | 伴随问题中的结果 |
| --- | --- | --- |
| $\mathbb C:\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$ | 用 $(\mathbb C:\boldsymbol A):\boldsymbol B=\boldsymbol A:(\mathbb C^T:\boldsymbol B)$ | $\mathbb C^T:\dot{\boldsymbol\varepsilon}_{\boldsymbol v}$ |
| 动量方程中的 $\nabla\tilde p$ | 分部积分后体内项为 $-\tilde p\,\nabla\cdot\boldsymbol v$ | $\nabla\cdot\boldsymbol v=0$ |
| 质量方程中的 $\nabla\cdot\tilde{\boldsymbol u}$ | 分部积分后体内项为 $\nabla r\cdot\tilde{\boldsymbol u}$ | 伴随应力中出现 $-r\boldsymbol I$ |
| 黏性体积分中的 $\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$ | 再分部积分，把导数转移到 $\mathbb C^T:\dot{\boldsymbol\varepsilon}_{\boldsymbol v}$ 上 | $-\nabla\cdot\boldsymbol\tau^*=0$ |
| 顶部自然牵引 $\boldsymbol\tau\boldsymbol n=0$ | 顶部伴随边界项与 $\delta J$ 相消 | $\boldsymbol\tau^*\boldsymbol n=-(u_x-u_{\rm obs})\boldsymbol e_x$ |
| 底部不可穿透 $\tilde{\boldsymbol u}\cdot\boldsymbol n=0$ | 只剩切向速度方向 | $\boldsymbol v\cdot\boldsymbol n=0$ |
| 底部滑移 $\beta\tilde{\boldsymbol u}_t$ | $\beta\boldsymbol T$ 自伴随 | $\beta\boldsymbol v_t$ |
| 周期条件 | 两侧边界项抵消 | 伴随速度和伴随牵引周期 |

## 10. 逻辑小结

连续伴随问题不是凭空写出来的。它来自下面的消元过程：

1. 用 $(\boldsymbol v,r)$ 与增量正问题弱配对；
2. 对体积分分部积分，把导数从 $(\tilde{\boldsymbol u},\tilde p)$ 转移到 $(\boldsymbol v,r)$；
3. 令任意的 $\tilde p$ 系数为零，得到 $\nabla\cdot\boldsymbol v=0$；
4. 令任意的 $\tilde{\boldsymbol u}$ 体积分系数为零，得到 $-\nabla\cdot\boldsymbol\tau^*=0$；
5. 用顶部边界项抵消目标函数变分，得到顶部伴随牵引；
6. 用底部边界项消去 $\tilde{\boldsymbol u}_t$，剩下 $\delta\beta$ 项，得到梯度公式。

因此，伴随方程可以理解为“增量正问题线性化算子的转置问题”，而梯度公式来自这个转置过程在底部滑移边界上留下的参数项。
