# 增量伴随问题的推导

本文解释 [equations.md](equations.md) 中增量伴随问题的来源。它用于 Gauss--Newton 方法中计算 Hessian-vector product：给定参数方向 $\delta\beta$，先求状态增量 $(\tilde{\boldsymbol u},\tilde p)$，再求伴随变量沿这个方向的一阶变化 $(\tilde{\boldsymbol v},\tilde r)$。

## 1. 从伴随问题开始

当前参数 $\beta$ 下，正问题解为 $(\boldsymbol u,p)$，伴随解为 $(\boldsymbol v,r)$。伴随问题可以抽象写成

$$
\mathcal A(\boldsymbol u,\beta)^*(\boldsymbol v,r)
=
-J_{\boldsymbol u}(\boldsymbol u),
$$

其中 $\mathcal A(\boldsymbol u,\beta)$ 是增量正问题的状态线性化算子，$\mathcal A(\boldsymbol u,\beta)^*$ 是它的伴随算子。右端 $J_{\boldsymbol u}(\boldsymbol u)$ 是目标函数对状态的一阶导数。

对本文的顶部水平速度观测，

$$
J(\boldsymbol u)
=
\frac12
\int_{\Gamma_t}
(u_x-u_{\rm obs})^2\,ds,
$$

所以

$$
J_{\boldsymbol u}(\boldsymbol u)[\boldsymbol w]
=
\int_{\Gamma_t}
(u_x-u_{\rm obs})w_x\,ds.
$$

这正对应伴随问题顶部边界条件

$$
\delta\boldsymbol\sigma^*(\boldsymbol u,p)[\boldsymbol v,r]\boldsymbol n
=
-(u_x-u_{\rm obs})\boldsymbol e_x
\qquad\text{on }\Gamma_t.
$$

## 2. 沿参数方向线性化

给定参数方向 $\delta\beta$，考虑一族扰动参数

$$
\beta_\epsilon
=
\beta+\epsilon\delta\beta.
$$

对应的正问题解和伴随解记为

$$
(\boldsymbol u_\epsilon,p_\epsilon)
=
(\boldsymbol u(\beta_\epsilon),p(\beta_\epsilon)),
\qquad
(\boldsymbol v_\epsilon,r_\epsilon)
=
(\boldsymbol v(\beta_\epsilon),r(\beta_\epsilon)).
$$

这里 $(\boldsymbol v_\epsilon,r_\epsilon)$ 不是任意函数，而是参数取
$\beta_\epsilon$ 时重新求得的伴随解。伴随问题本身依赖当前参数和当前正问题解：

$$
\beta
\longmapsto
(\boldsymbol u(\beta),p(\beta))
\longmapsto
\mathcal A(\boldsymbol u(\beta),\beta)^*,
\qquad
\beta
\longmapsto
\boldsymbol u(\beta)
\longmapsto
J_{\boldsymbol u}(\boldsymbol u(\beta)).
$$

因此，当参数从 $\beta$ 变成 $\beta_\epsilon$ 时，伴随方程也要在新的状态
$(\boldsymbol u_\epsilon,p_\epsilon)$ 和新的参数 $\beta_\epsilon$ 处重新写出。

它们对 $\epsilon$ 的导数定义为

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

对每个 $\epsilon$，伴随方程都满足

$$
\mathcal A(\boldsymbol u_\epsilon,\beta_\epsilon)^*
(\boldsymbol v_\epsilon,r_\epsilon)
=
-J_{\boldsymbol u}(\boldsymbol u_\epsilon).
$$

也就是说，

$$
F(\epsilon)
:=
\mathcal A(\boldsymbol u_\epsilon,\beta_\epsilon)^*
(\boldsymbol v_\epsilon,r_\epsilon)
+
J_{\boldsymbol u}(\boldsymbol u_\epsilon)
=0
$$

对每个 $\epsilon$ 都成立。增量伴随方程就是对这个恒等式在
$\epsilon=0$ 处求导得到的。

现在对这个等式求导。左端是“算子作用在变量上”的复合函数，所以有两部分：

$$
\left.
\frac{d}{d\epsilon}
\left[
\mathcal A(\boldsymbol u_\epsilon,\beta_\epsilon)^*
(\boldsymbol v_\epsilon,r_\epsilon)
\right]
\right|_{\epsilon=0}
=
\mathcal A(\boldsymbol u,\beta)^*
(\tilde{\boldsymbol v},\tilde r)
+
\delta\!\left(\mathcal A(\boldsymbol u,\beta)^*\right)
[\tilde{\boldsymbol u},\delta\beta](\boldsymbol v,r).
$$

第一项来自伴随未知量本身的变化：

$$
(\boldsymbol v,r)
\longmapsto
(\boldsymbol v,r)
+
\epsilon(\tilde{\boldsymbol v},\tilde r).
$$

第二项来自伴随算子系数的变化。因为 $\mathcal A^*$ 依赖当前正问题解和参数，即依赖 $(\boldsymbol u,p,\beta)$，所以当

$$
(\boldsymbol u,p,\beta)
\longmapsto
(\boldsymbol u,p,\beta)
+
\epsilon(\tilde{\boldsymbol u},\tilde p,\delta\beta)
$$

时，算子本身也会改变。这部分记为

$$
\delta\!\left(\mathcal A(\boldsymbol u,\beta)^*\right)
[\tilde{\boldsymbol u},\delta\beta](\boldsymbol v,r).
$$

它表示“先把伴随算子沿 $(\tilde{\boldsymbol u},\delta\beta)$ 方向线性化，再作用到当前伴随解 $(\boldsymbol v,r)$ 上”。

右端只通过状态 $\boldsymbol u_\epsilon$ 变化：

$$
\left.
\frac{d}{d\epsilon}
\left[
-J_{\boldsymbol u}(\boldsymbol u_\epsilon)
\right]
\right|_{\epsilon=0}
=
-J_{\boldsymbol u\boldsymbol u}(\boldsymbol u)
[\tilde{\boldsymbol u},\cdot].
$$

这里的 $\cdot$ 表示这个二阶变分仍然是一个线性泛函，要作用在测试速度或伴随弱式中的边界速度方向上。

把左端和右端合并，得到

$$
\mathcal A(\boldsymbol u,\beta)^*
(\tilde{\boldsymbol v},\tilde r)
+
\delta\!\left(\mathcal A(\boldsymbol u,\beta)^*\right)
[\tilde{\boldsymbol u},\delta\beta](\boldsymbol v,r)
=
-J_{\boldsymbol u\boldsymbol u}(\boldsymbol u)[\tilde{\boldsymbol u},\cdot].
$$

这就是完整增量伴随方程的抽象形式。第一项是用原来的伴随算子作用在增量伴随变量上；第二项是伴随算子本身随 $(\boldsymbol u,\beta)$ 改变而产生的变化；右端是目标函数二阶变分。

## 3. Gauss--Newton 近似

`equations.md` 中写的是 Gauss--Newton 增量伴随问题。Gauss--Newton 近似保留观测项的二阶变化

$$
J_{\boldsymbol u\boldsymbol u}(\boldsymbol u)[\tilde{\boldsymbol u},\cdot],
$$

但忽略伴随算子本身的变化项

$$
\delta\!\left(\mathcal A(\boldsymbol u,\beta)^*\right)
[\tilde{\boldsymbol u},\delta\beta](\boldsymbol v,r).
$$

因此近似方程变成

$$
\mathcal A(\boldsymbol u,\beta)^*
(\tilde{\boldsymbol v},\tilde r)
=
-J_{\boldsymbol u\boldsymbol u}(\boldsymbol u)[\tilde{\boldsymbol u},\cdot].
$$

这说明增量伴随问题使用和伴随问题相同的左端算子，只是右端从目标函数一阶变分换成目标函数二阶变分。

## 4. 顶部源项为什么是 $-\tilde u_x\boldsymbol e_x$

目标函数一阶变分是

$$
J_{\boldsymbol u}(\boldsymbol u)[\boldsymbol w]
=
\int_{\Gamma_t}
(u_x-u_{\rm obs})w_x\,ds.
$$

沿状态方向 $\tilde{\boldsymbol u}$ 再求一次变分：

$$
J_{\boldsymbol u\boldsymbol u}(\boldsymbol u)
[\tilde{\boldsymbol u},\boldsymbol w]
=
\int_{\Gamma_t}
\tilde u_x w_x\,ds.
$$

也就是

$$
J_{\boldsymbol u\boldsymbol u}(\boldsymbol u)
[\tilde{\boldsymbol u},\boldsymbol w]
=
\int_{\Gamma_t}
\tilde u_x\boldsymbol e_x
\cdot\boldsymbol w\,ds.
$$

由于伴随方程采用的符号约定是顶部伴随牵引等于负的目标函数状态导数，所以增量伴随顶部牵引为

$$
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
=
-\tilde u_x\boldsymbol e_x
\qquad\text{on }\Gamma_t.
$$

这就是 `equations.md` 中增量伴随顶部边界条件的来源。

## 5. 强形式

在 Gauss--Newton 近似下，左端算子冻结在当前正问题解 $(\boldsymbol u,p)$ 处。因此增量伴随应力仍写成

$$
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]
=
\mathbb C(\boldsymbol u)^T:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol v}}
-
\tilde r\boldsymbol I.
$$

于是增量伴随问题为求 $(\tilde{\boldsymbol v},\tilde r)$，满足

$$
\left\{
\begin{aligned}
\nabla\cdot\tilde{\boldsymbol v} &=0
&&\text{in }\Omega,\\
-\nabla\cdot
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]
&=\boldsymbol 0
&&\text{in }\Omega,\\
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
&=
-\tilde u_x\boldsymbol e_x
&&\text{on }\Gamma_t,\\
\tilde{\boldsymbol v}\cdot\boldsymbol n &=0
&&\text{on }\Gamma_b,\\
\boldsymbol T
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
+
\beta\tilde{\boldsymbol v}_t
&=\boldsymbol 0
&&\text{on }\Gamma_b,\\
\tilde{\boldsymbol v},\
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
&\text{ periodic}
&&\text{on }\Gamma_p.
\end{aligned}
\right.
$$

它和伴随问题的区别只有顶部源项：

$$
-(u_x-u_{\rm obs})\boldsymbol e_x
\quad\longrightarrow\quad
-\tilde u_x\boldsymbol e_x.
$$

前者由当前观测误差驱动，后者由增量状态在观测量上的变化驱动。

## 6. 与完整 Hessian 的区别

如果不做 Gauss--Newton 近似，线性化伴随问题还会包含

$$
\delta\!\left(\mathcal A(\boldsymbol u,\beta)^*\right)
[\tilde{\boldsymbol u},\delta\beta](\boldsymbol v,r),
$$

也就是非线性 Stokes 切线算子、底部滑移系数和当前伴随解共同产生的附加项。例如底部伴随条件

$$
\boldsymbol T
\delta\boldsymbol\sigma^*(\boldsymbol u,p)[\boldsymbol v,r]\boldsymbol n
+
\beta\boldsymbol v_t
=\boldsymbol 0
$$

完整线性化时会产生

$$
\boldsymbol T
\delta\boldsymbol\sigma^*(\boldsymbol u,p)
[\tilde{\boldsymbol v},\tilde r]\boldsymbol n
+
\beta\tilde{\boldsymbol v}_t
$$

之外的项，例如由 $\delta\beta\,\boldsymbol v_t$ 和应力算子对 $\boldsymbol u$ 的变化引起的项。`equations.md` 中的增量伴随问题把这些项作为完整 Hessian 的非 Gauss--Newton 部分忽略掉，所以底部条件保持齐次。

因此，`equations.md` 中的增量伴随问题应理解为：固定当前线性化 Stokes 算子，只把顶部观测误差的变化传回到底部，用于构造 Gauss--Newton Hessian-vector product。
