# 通用滑移指数 m 的线性化

[equations.md](equations.md) 已整理线性滑移情形 `m=1` 下的正问题、增量正问题、伴随问题和增量伴随问题。本文只记录通用 Weertman 滑移指数 `m` 时需要替换的底部项，避免重复完整方程。

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

记速度扰动的底部切向分量为

$$
\tilde{\boldsymbol u}_t=\boldsymbol T\tilde{\boldsymbol u}.
$$

对底部摩擦项求速度方向导数：

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

`m=1` 时第二项为零，第一项退化为

$$
\beta\tilde{\boldsymbol u}_t,
$$

这正是 [equations.md](equations.md) 中使用的线性滑移形式。

## 对参数的线性化

对参数扰动 `\delta\beta`，有

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
D_{\boldsymbol u}
\left(
\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t
\right)
[\tilde{\boldsymbol u}],
$$

并把右端参数项

$$
-\delta\beta\,\boldsymbol u_t
$$

替换为

$$
-\delta\beta|\boldsymbol u_t|^{m-1}\boldsymbol u_t.
$$

## 伴随问题中的替换

伴随底部切向条件中使用速度线性化算子的伴随。由于

$$
\boldsymbol u_t\mapsto
|\boldsymbol u_t|^{m-1}\boldsymbol u_t
$$

的 Jacobian 为

$$
|\boldsymbol u_t|^{m-1}\boldsymbol I_t
+
(m-1)|\boldsymbol u_t|^{m-3}
\boldsymbol u_t\otimes\boldsymbol u_t,
$$

其中 `\boldsymbol I_t` 表示切向空间上的恒等映射。该算子在通常欧氏内积下是对称的，所以伴随底部项形式上仍为同一个切向线性化算子作用在伴随速度切向分量上：

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

`m=1` 时该项退化为

$$
\beta\boldsymbol v_t.
$$

## 梯度公式

通用 `m` 时，目标函数对 `\beta` 的连续梯度密度为

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

## 数值注意

当 `m<1` 或底部速度接近零时，`|\boldsymbol u_t|^{m-3}` 可能带来奇异或病态系数。实现中通常需要与正问题相同的正则化，例如把

$$
|\boldsymbol u_t|
$$

替换为

$$
\left(|\boldsymbol u_t|^2+\varepsilon_{\mathrm{reg}}^2\right)^{1/2}.
$$

离散装配时，正问题残差、增量正问题、伴随问题和梯度必须使用同一套正则化，否则伴随梯度会和有限差分检查不一致。
