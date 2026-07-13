# 非线性 Stokes 的严格残差、乘子阻尼与梯度检查

本文解释当前非线性 Stokes 正问题和反问题实现中的三个问题：

1. 什么是严格残差检查；
2. 什么是拉格朗日乘子与速度、压力同步阻尼；
3. 当前是否验证了伴随梯度与有限差分梯度的一致性。

相关代码：

```text
ice-sheet/FullStokes2d/NonlinearStokesP2P1.m
ice-sheet/FullStokes2d/NSAdjInvTikhonov.m
ice-sheet/FullStokes2d/NSFDInversion.m
ice-sheet/FullStokes2d/unused/NSRegression.m
```

## 1. 非线性离散方程

非线性 Stokes 有限元方程可以写成：

$$
R(u,p,\lambda)=0,
$$

其中：

- $u$ 是速度自由度；
- $p$ 是压力自由度；
- $\lambda$ 是用于周期条件、底部不可穿透条件和压力规范的拉格朗日乘子。

把速度和压力组合为

$$
x=
\begin{bmatrix}
u\\
p
\end{bmatrix},
$$

离散状态矩阵为

$$
M(u)=
\begin{bmatrix}
K(u)+K_b(u) & B^T\\
B & 0
\end{bmatrix}.
$$

这里：

- $K(u)$ 是 Glen 非线性黏性矩阵；
- $K_b(u)$ 是非线性底部滑移矩阵；
- $B$ 是离散散度矩阵；
- $C$ 是周期条件、底部不可穿透条件和压力规范的约束矩阵。

完整增广系统为

$$
\begin{bmatrix}
M(u) & C^T\\
C & 0
\end{bmatrix}
\begin{bmatrix}
x\\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
b\\
0
\end{bmatrix},
$$

其中

$$
b=
\begin{bmatrix}
F\\
0
\end{bmatrix}.
$$

## 2. 只检查迭代增量有什么不足

原来的 Picard 停止条件主要检查相邻两次速度的相对变化：

$$
r_{\mathrm{change}}^k
=
\frac{\|u^{k+1}-u^k\|}
{\max(1,\|u^{k+1}\|)}.
$$

当

$$
r_{\mathrm{change}}^k<\mathrm{tol}
$$

时，就认为迭代收敛。

这个条件只能说明：

> 当前迭代已经不再明显改变。

它不能单独证明：

> 当前状态真正满足非线性离散方程。

以下情况可能导致增量很小，但方程残差仍然不小。

### 2.1 阻尼很小

如果阻尼参数 $\alpha$ 很小，

$$
u^{k+1}=(1-\alpha)u^k+\alpha\widehat u^{k+1},
$$

那么即使未阻尼更新 $\widehat u^{k+1}-u^k$ 很大，实际更新

$$
u^{k+1}-u^k
=\alpha(\widehat u^{k+1}-u^k)
$$

仍可能很小。

### 2.2 Picard 迭代停滞

迭代可能进入一个变化缓慢但没有真正满足方程的状态。此时只看相邻解差不能识别停滞。

### 2.3 参数反演中的热启动

参数从 $\beta^k$ 变为 $\beta^{k+1}$ 后，程序通常使用上一次正问题解作为初值。

该初值可能与新参数对应的解很接近，因此第一步速度变化很小，但它不一定已经满足新参数下的方程：

$$
R(u^k,p^k,\lambda^k;\beta^{k+1})\ne0.
$$

如果只检查速度增量，可能过早停止。

## 3. 什么是严格残差检查

“严格残差检查”是指在更新后的当前状态重新装配非线性算子，并直接检查完整离散方程。

在当前速度 $u$ 上重新装配：

$$
K(u),\qquad K_b(u),\qquad F.
$$

状态残差定义为

$$
R_{\mathrm{state}}
=M(u)x-b+C^T\lambda.
$$

约束残差定义为

$$
R_{\mathrm{constraint}}=Cx.
$$

完整增广残差为

$$
R_{\mathrm{aug}}
=
\begin{bmatrix}
R_{\mathrm{state}}\\
R_{\mathrm{constraint}}
\end{bmatrix}.
$$

归一化状态残差为

$$
r_{\mathrm{state}}
=
\frac{\|R_{\mathrm{state}}\|}
{\max(1,\|b\|)}.
$$

归一化约束残差为

$$
r_{\mathrm{constraint}}
=
\frac{\|Cx\|}
{\max(1,\|x\|)}.
$$

程序使用

$$
r_{\mathrm{nonlinear}}
=
\max\left(
r_{\mathrm{state}},
r_{\mathrm{constraint}}
\right)
$$

作为完整非线性残差。

现在只有同时满足

$$
r_{\mathrm{change}}<\mathrm{tol}
$$

和

$$
r_{\mathrm{nonlinear}}
<\mathrm{residual\_tol}
$$

时，求解器才会设置

```matlab
info.converged = true;
```

## 4. 为什么称为“完整”残差

这个残差包含：

- 当前速度下重新计算的 Glen 黏度；
- 当前速度下重新计算的底部滑移系数；
- 压力梯度；
- 不可压缩方程；
- 体力和边界牵引；
- 周期约束反力；
- 底部不可穿透约束反力；
- 如果启用，压力规范对应的约束反力。

它检查的是当前有限元离散系统，而不是只检查某个冻结系数的 Picard 线性系统。

需要注意，它仍然是离散残差，不是连续 PDE 在每个空间点的强形式残差。

## 5. 为什么需要拉格朗日乘子

约束写成

$$
Cx=0.
$$

如果只解状态方程而不引入乘子，无法在保持对称鞍点结构的同时强制这些约束。

引入拉格朗日乘子后，状态方程中出现约束反力：

$$
M(u)x+C^T\lambda=b.
$$

$\lambda$ 可以理解为维持约束所需的离散反力。例如：

- 周期边界对应的匹配反力；
- 底部不可穿透条件对应的法向反力；
- 压力零均值规范对应的乘子。

因此计算状态残差时必须保留

$$
C^T\lambda.
$$

如果漏掉这一项，即使增广方程已经精确满足，也会把合法的约束反力误判为动量或压力方程误差。

## 6. 什么是同步阻尼

在一次 Picard 线性求解中，程序得到未阻尼候选解：

$$
\widehat u^{k+1},
\qquad
\widehat p^{k+1},
\qquad
\widehat\lambda^{k+1}.
$$

为了提高非线性迭代稳定性，速度和压力不直接替换，而是使用阻尼：

$$
u^{k+1}
=(1-\alpha)u^k+\alpha\widehat u^{k+1},
$$

$$
p^{k+1}
=(1-\alpha)p^k+\alpha\widehat p^{k+1}.
$$

其中

$$
0<\alpha\le1.
$$

同步阻尼表示拉格朗日乘子采用同样的更新：

$$
\lambda^{k+1}
=(1-\alpha)\lambda^k
+\alpha\widehat\lambda^{k+1}.
$$

所以当前迭代状态是统一的三元组：

$$
\left(
u^{k+1},
p^{k+1},
\lambda^{k+1}
\right).
$$

## 7. 为什么乘子必须同步阻尼

完整状态残差包含

$$
M(u)x+C^T\lambda-b.
$$

如果速度和压力使用阻尼后的值，但乘子直接使用未阻尼候选值：

$$
x^{k+1}
=(1-\alpha)x^k+\alpha\widehat x^{k+1},
$$

$$
\lambda^{k+1}
=\widehat\lambda^{k+1},
$$

那么 $x^{k+1}$ 和 $\lambda^{k+1}$ 来自两个不同的迭代状态。

这会使残差混入人为的不一致：

$$
M(u^{k+1})x^{k+1}
+C^T\widehat\lambda^{k+1}-b.
$$

即使迭代本身正在正常收敛，该量也不能正确反映阻尼后状态的方程误差。

同步阻尼后，速度、压力和约束反力沿同一个步长更新，残差检查才具有一致含义。

## 8. 严格残差检查带来的计算成本

严格残差需要在更新后的速度上重新装配：

$$
K(u^{k+1}),\qquad K_b(u^{k+1}).
$$

这会增加计算量。

当前实现采用延迟检查：

1. 先计算每步迭代增量；
2. 当增量低于 `residual_check_threshold` 时开始完整残差检查；
3. 一旦开始检查，之后每步都检查；
4. 最大迭代步的最后一步始终检查。

这样可以避免在离收敛还很远时重复装配完整残差。

有限差分反演中，热启动正问题原先可能只用一步就因速度增量很小而退出。加入完整残差后，通常仍需约 15 步，直到新参数下的离散方程真正满足容差。

## 9. 当前伴随梯度验证了什么

设反演变量为

$$
q=\log\beta,
$$

目标函数为

$$
J(q).
$$

伴随法计算梯度：

$$
g_{\mathrm{adj}}=\nabla_qJ.
$$

程序选取一个归一化方向 $d$，通过中心有限差分计算方向导数：

$$
D_{\mathrm{FD}}
=
\frac{
J(q+h d)-J(q-h d)
}{
2h
}.
$$

同时计算伴随梯度方向导数：

$$
D_{\mathrm{adj}}
=g_{\mathrm{adj}}^Td.
$$

然后比较相对误差：

$$
e_g
=
\frac{
|D_{\mathrm{FD}}-D_{\mathrm{adj}}|
}{
\max\left(
\epsilon,
|D_{\mathrm{FD}}|,
|D_{\mathrm{adj}}|
\right)
}.
$$

最近一次完整回归测试得到：

```text
finite-difference direction = 5.274947e-2
adjoint direction           = 5.274938e-2
relative gradient error     = 1.551e-6
```

因此，在同一个参数初值和同一个参数方向上，伴随梯度与目标函数中心有限差分方向导数一致。

## 10. 同时还验证了哪些导数

伴随脚本还检查了增量状态。

一致切线方程给出：

$$
R_U\delta U=-R_qd.
$$

有限差分状态方向为

$$
\delta U_{\mathrm{FD}}
=
\frac{
U(q+hd)-U(q-hd)
}{
2h
}.
$$

程序比较：

$$
\delta U
\quad\text{和}\quad
\delta U_{\mathrm{FD}}.
$$

最近结果为：

```text
state derivative error = 2.396e-7
```

程序还检查了 Gauss--Newton 方向量，最近结果为：

```text
Gauss-Newton check error = 4.797e-7
```

这些结果说明当前切线矩阵、参数导数和伴随方程在选定方向上相互一致。

## 11. 当前尚未完成的梯度比较

当前已经完成的是：

$$
g_{\mathrm{adj}}^Td
\quad\text{与}\quad
\frac{J(q+hd)-J(q-hd)}{2h}
$$

在一个方向 $d$ 上的比较。

尚未完成的是整个梯度向量逐分量比较：

$$
\left(g_{\mathrm{FD}}\right)_j
=
\frac{
J(q+h e_j)-J(q-h e_j)
}{
2h
},
$$

以及

$$
\frac{
\|g_{\mathrm{adj}}-g_{\mathrm{FD}}\|
}{
\max(
\|g_{\mathrm{adj}}\|,
\|g_{\mathrm{FD}}\|
)
}.
$$

另外，现有两个反演脚本不能直接把各自内部的梯度拿来横向比较，因为它们的设置并不相同，包括：

- 几何高度不同；
- Weertman 指数不同；
- 参数真值和初值不同；
- 观测分量不同；
- 正则化设置不同。

所以，“两个脚本都能运行”不等于“它们计算的是同一个目标函数梯度”。

## 12. 如何进行完整梯度向量检查

应建立一个统一测试问题，固定：

- 同一网格和几何；
- 同一个 Glen 指数；
- 同一个 Weertman 指数；
- 同一个参数初值 $q$；
- 同一组观测自由度；
- 同一个数据尺度；
- 同一个目标函数；
- 同一个正则化项；
- 同一个正问题容差。

然后按以下步骤检查。

### 12.1 计算伴随梯度

求正问题：

$$
R(U,q)=0.
$$

求伴随问题：

$$
R_U^T\psi=-J_U^T.
$$

计算：

$$
g_{\mathrm{adj}}
=R_q^T\psi+J_q.
$$

### 12.2 逐参数中心差分

对于每个单位方向 $e_j$：

$$
\left(g_{\mathrm{FD}}\right)_j
=
\frac{
J(q+h e_j)-J(q-h e_j)
}{
2h
}.
$$

### 12.3 比较完整向量

报告：

$$
e_{\mathrm{vector}}
=
\frac{
\|g_{\mathrm{adj}}-g_{\mathrm{FD}}\|
}{
\max(
\epsilon,
\|g_{\mathrm{adj}}\|,
\|g_{\mathrm{FD}}\|
)
}.
$$

同时报告每一分量：

| 参数分量 | 伴随梯度 | 有限差分梯度 | 相对差异 |
|---|---:|---:|---:|
| 1 |  |  |  |
| 2 |  |  |  |
| $\cdots$ |  |  |  |

这比单方向检查更容易定位某个参数基函数、周期插值或底部积分中的局部错误。

## 13. 当前结论

### 严格残差

当前求解器不再只依赖相邻 Picard 解的变化，而是同时要求完整增广非线性离散残差达到容差。

### 同步阻尼

速度、压力和约束拉格朗日乘子使用同一个阻尼参数更新，从而使残差中的状态和约束反力属于同一个迭代点。

### 梯度一致性

当前已经验证：

- 增量状态与中心有限差分一致；
- 伴随梯度方向导数与目标函数中心有限差分一致；
- Gauss--Newton 方向量与有限差分一致。

当前尚未验证：

- 同一统一测试问题下，伴随梯度向量与逐参数中心有限差分梯度向量的完整逐分量一致性。
