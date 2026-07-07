# 通用滑移指数 m 的线性化

[equations.md](equations.md) 已整理线性滑移情形 `m=1` 下的正问题、增量正问题、伴随问题和增量伴随问题。本文只记录通用 Weertman 滑移指数 `m` 时需要替换的底部项，避免重复完整方程。

本文的目的不是改变体内非线性 Stokes 方程，而是说明当底部滑移律从线性形式

$$
\beta\boldsymbol u_t
$$

换成幂律形式

$$
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
$$

以后，增量正问题、伴随问题和梯度公式中哪些项随之改变。

下面所有线性化都在当前正问题解 $(\boldsymbol u,p)$ 和当前参数
$\beta$ 处进行。速度方向记为 $\tilde{\boldsymbol u}$，参数方向记为
$\delta\beta$。

## 底部滑移律

线性滑移 `m=1` 时，底部摩擦项为

$$
\beta \boldsymbol u_t.
$$

通用 `m` 时替换为

$$
\beta |\boldsymbol u_t|^{m-1}\boldsymbol u_t,
$$

其中

$$
\boldsymbol u_t=\boldsymbol T\boldsymbol u,
\qquad
\boldsymbol T=\boldsymbol I-\boldsymbol n\otimes\boldsymbol n.
$$

因此底部切向边界条件为

$$
\boldsymbol T\boldsymbol\sigma\boldsymbol n
+
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
=\boldsymbol 0
\qquad\text{on }\Gamma_b.
$$

## 对速度的线性化

把底部摩擦项写成

$$
\boldsymbol b(\boldsymbol u,\beta)
=
\beta\boldsymbol \phi(\boldsymbol u_t),
\qquad
\boldsymbol \phi(\boldsymbol z)
=
|\boldsymbol z|^{m-1}\boldsymbol z.
$$

这里 $\boldsymbol z$ 是切向向量。因为 $\boldsymbol u_t=\boldsymbol T\boldsymbol u$，所以速度扰动对应的切向扰动是

$$
\tilde{\boldsymbol u}_t=\boldsymbol T\tilde{\boldsymbol u}.
$$

先求 $\boldsymbol \phi$ 对 $\boldsymbol z$ 的方向导数。令

$$
\boldsymbol z_\epsilon
=
\boldsymbol z+\epsilon\boldsymbol h.
$$

则

$$
\left.
\frac{d}{d\epsilon}
\left(
|\boldsymbol z_\epsilon|^{m-1}
\boldsymbol z_\epsilon
\right)
\right|_{\epsilon=0}
=
|\boldsymbol z|^{m-1}\boldsymbol h
+
\left.
\frac{d}{d\epsilon}
|\boldsymbol z_\epsilon|^{m-1}
\right|_{\epsilon=0}
\boldsymbol z.
$$

又因为

$$
\left.
\frac{d}{d\epsilon}
|\boldsymbol z+\epsilon\boldsymbol h|^2
\right|_{\epsilon=0}
=
2\boldsymbol z\cdot\boldsymbol h,
$$

所以

$$
\left.
\frac{d}{d\epsilon}
|\boldsymbol z+\epsilon\boldsymbol h|^{m-1}
\right|_{\epsilon=0}
=
(m-1)|\boldsymbol z|^{m-3}
\boldsymbol z\cdot\boldsymbol h.
$$

因此

$$
D\boldsymbol \phi(\boldsymbol z)[\boldsymbol h]
=
|\boldsymbol z|^{m-1}\boldsymbol h
+
(m-1)|\boldsymbol z|^{m-3}
(\boldsymbol z\cdot\boldsymbol h)\boldsymbol z.
$$

令 $\boldsymbol z=\boldsymbol u_t$、$\boldsymbol h=\tilde{\boldsymbol u}_t$，并乘上固定参数 $\beta$，得到底部摩擦项对速度的方向导数：

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
&\quad+
\beta(m-1)|\boldsymbol u_t|^{m-3}
\left(\boldsymbol u_t\cdot\tilde{\boldsymbol u}_t\right)
\boldsymbol u_t.
\end{aligned}
$$

这个表达式可写成一个切向线性算子

$$
\boldsymbol K_b(\boldsymbol u,\beta)\tilde{\boldsymbol u}_t,
$$

其中

$$
\boldsymbol K_b(\boldsymbol u,\beta)
=
\beta
\left[
|\boldsymbol u_t|^{m-1}\boldsymbol I_t
+
(m-1)|\boldsymbol u_t|^{m-3}
\boldsymbol u_t\otimes\boldsymbol u_t
\right].
$$

$\boldsymbol I_t$ 表示切向空间上的恒等映射。由于
$\boldsymbol u_t\otimes\boldsymbol u_t$ 是对称二阶张量，$\boldsymbol K_b$ 在切向欧氏内积下是自伴随的：

$$
\boldsymbol K_b^*=\boldsymbol K_b.
$$

`m=1` 时第二项为零，第一项退化为

$$
\boldsymbol K_b\tilde{\boldsymbol u}_t
=
\beta\tilde{\boldsymbol u}_t,
$$

这正是 [equations.md](equations.md) 中使用的线性滑移形式。

## 对参数的线性化

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

也就是说，参数方向只给出一个底部已知源项；它不改变体内动量方程和不可压缩方程。

如果使用对数参数

$$
q=\log\beta,
$$

则

$$
\delta\beta=\beta\,\delta q.
$$

因此参数方向导数可写成

$$
\beta\,\delta q
|\boldsymbol u_t|^{m-1}\boldsymbol u_t.
$$

## 增量正问题中的替换

在 [equations.md](equations.md) 的增量正问题底部切向条件中，把线性滑移项

$$
\beta\tilde{\boldsymbol u}_t
$$

替换为

$$
\boldsymbol K_b(\boldsymbol u,\beta)\tilde{\boldsymbol u}_t,
$$

并把右端参数项

$$
-\delta\beta\,\boldsymbol u_t
$$

替换为

$$
-\delta\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t.
$$

因此通用 `m` 的底部增量条件为

$$
\boldsymbol T
\delta\boldsymbol\sigma(\boldsymbol u,p)
[\tilde{\boldsymbol u},\tilde p]\boldsymbol n
+
\boldsymbol K_b(\boldsymbol u,\beta)\tilde{\boldsymbol u}_t
=
-\delta\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\qquad\text{on }\Gamma_b.
$$

其余边界条件和体内方程与 [equations.md](equations.md) 相同。

## 伴随问题中的替换

伴随底部切向条件中使用速度线性化算子的伴随。上面已经得到

$$
\boldsymbol K_b^*=\boldsymbol K_b,
$$

所以伴随底部项形式上仍为同一个切向线性化算子作用在伴随速度切向分量上：

$$
\boldsymbol K_b^*\boldsymbol v_t
=
\boldsymbol K_b\boldsymbol v_t.
$$

`m=1` 时该项退化为

$$
\beta\boldsymbol v_t.
$$

因此通用 `m` 的伴随底部切向条件为

$$
\boldsymbol T
\delta\boldsymbol\sigma^*(\boldsymbol u,p)[\boldsymbol v,r]\boldsymbol n
+
\boldsymbol K_b(\boldsymbol u,\beta)\boldsymbol v_t
=
\boldsymbol 0
\qquad\text{on }\Gamma_b.
$$

增量伴随问题在 Gauss--Newton 近似下使用相同冻结算子，因此只需把
$\boldsymbol v_t$ 换成 $\tilde{\boldsymbol v}_t$。

## 梯度公式

伴随消元以后，底部只剩由参数扰动产生的项。由于参数残差为

$$
\delta\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t,
$$

它与伴随切向速度配对得到

$$
\delta J(\beta)[\delta\beta]
=
\int_{\Gamma_b}
\delta\beta
|\boldsymbol u_t|^{m-1}
\boldsymbol u_t\cdot\boldsymbol v_t
\,ds.
$$

因此通用 `m` 时，目标函数对 `\beta` 的连续梯度密度为

$$
g_\beta
=
|\boldsymbol u_t|^{m-1}
\boldsymbol u_t\cdot\boldsymbol v_t.
$$

对数参数 `q=\log\beta` 下，

$$
g_q
=
\beta
|\boldsymbol u_t|^{m-1}
\boldsymbol u_t\cdot\boldsymbol v_t.
$$

`m=1` 时退化为

$$
g_q=\beta\boldsymbol u_t\cdot\boldsymbol v_t.
$$

符号仍依赖伴随方程中顶部牵引的约定。本文沿用
[equations.md](equations.md) 和 [adjoint.md](adjoint.md) 的约定，因此梯度符号与那两篇一致。

## 数值注意

当 `m<1` 或底部速度接近零时，$|\boldsymbol u_t|^{m-3}$ 可能带来奇异或病态系数。实现中通常需要与正问题相同的正则化，例如把

$$
|\boldsymbol u_t|
$$

替换为

$$
\left(|\boldsymbol u_t|^2+\varepsilon_{\mathrm{reg}}^2\right)^{1/2}.
$$

离散装配时，正问题残差、增量正问题、伴随问题和梯度必须使用同一套正则化，否则伴随梯度会和有限差分检查不一致。

若采用正则化

$$
s_\varepsilon
=
\left(|\boldsymbol u_t|^2+\varepsilon_{\mathrm{reg}}^2\right)^{1/2},
$$

则滑移律通常写成

$$
\beta s_\varepsilon^{m-1}\boldsymbol u_t.
$$

此时速度线性化中的 Jacobian 应同步改为

$$
\beta
\left[
s_\varepsilon^{m-1}\boldsymbol I_t
+
(m-1)s_\varepsilon^{m-3}
\boldsymbol u_t\otimes\boldsymbol u_t
\right].
$$

注意这里第二项仍是 $\boldsymbol u_t\otimes\boldsymbol u_t$，不是
$s_\varepsilon^2\boldsymbol I_t$。正则化只改变幂函数的分母尺度，不改变
$|\boldsymbol u_t|^2$ 对速度方向的导数结构。
