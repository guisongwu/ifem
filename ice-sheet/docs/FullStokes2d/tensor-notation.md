# 非线性 Stokes 推导中的张量记号

本文单独整理 `equations.md` 中用到的张量记号，包括张量积、二阶张量双点积、四阶张量、四阶张量作用在二阶张量上，以及四阶张量的转置。

## 张量积

对任意向量 $\boldsymbol a,\boldsymbol b\in\mathbb R^d$，张量积 $\boldsymbol a\otimes\boldsymbol b$ 是一个二阶张量，其分量定义为

$$
(\boldsymbol a\otimes\boldsymbol b)_{ij}
=
a_i b_j.
$$

在底部边界上，若 $\boldsymbol n$ 是单位外法向，则

$$
\boldsymbol n\otimes\boldsymbol n
$$

给出法向投影。对应的切向投影算子为

$$
\boldsymbol T
=
\boldsymbol I-\boldsymbol n\otimes\boldsymbol n.
$$

因此任意速度 $\boldsymbol w$ 的切向分量可以写成

$$
\boldsymbol w_t
=
\boldsymbol T\boldsymbol w.
$$

## 二阶张量的双点积

两个二阶张量 $\boldsymbol A,\boldsymbol B$ 的双点积定义为

$$
\boldsymbol A:\boldsymbol B
=
A_{ij}B_{ij}.
$$

这里默认对重复指标 $i,j$ 求和。二维时就是

$$
\boldsymbol A:\boldsymbol B
=
A_{11}B_{11}
+
A_{12}B_{12}
+
A_{21}B_{21}
+
A_{22}B_{22}.
$$

在 Stokes 推导中，应变率第二不变量写成

$$
\varepsilon_{\rm II}
=
\frac12
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}.
$$

如果速度沿方向 $\tilde{\boldsymbol u}$ 扰动，则应变率扰动为 $\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$，于是

$$
\delta\varepsilon_{\rm II}(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}.
$$

因此当黏度只依赖 $\varepsilon_{\rm II}$ 时，由链式法则有

$$
\eta'(\boldsymbol u)[\tilde{\boldsymbol u}]
=
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right).
$$

## 四阶张量

四阶张量可以理解为带四个指标的量：

$$
C_{ijkl}.
$$

它在这里的作用是表示“二阶张量到二阶张量”的线性映射。类比矩阵：

$$
y_i=A_{ij}x_j
$$

表示矩阵 $A$ 把向量 $x$ 映射成向量 $y$。如果输入和输出都换成二阶张量，那么映射就写成

$$
B_{ij}=C_{ijkl}A_{kl}.
$$

这里 $C_{ijkl}$ 就是四阶张量。

## 四阶张量与二阶张量的双点积

若 $\mathbb C$ 是四阶张量，$\boldsymbol A$ 是二阶张量，则 $\mathbb C:\boldsymbol A$ 是一个二阶张量，定义为

$$
\left(\mathbb C:\boldsymbol A\right)_{ij}
=
C_{ijkl}A_{kl}.
$$

也就是说，对重复指标 $k,l$ 求和，剩下的 $i,j$ 组成结果二阶张量。

二维中，例如

$$
\left(\mathbb C:\boldsymbol A\right)_{11}
=
C_{1111}A_{11}
+
C_{1112}A_{12}
+
C_{1121}A_{21}
+
C_{1122}A_{22}.
$$

这和矩阵乘向量完全类比，只是矩阵的一个输入指标变成了两个输入指标。

## 为什么黏性切线可以写成四阶张量

非线性黏性应力的黏性部分为

$$
\boldsymbol\sigma_{\rm visc}
=
2\eta(\boldsymbol u)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}.
$$

线性化以后，黏性应力的一阶变化为

$$
\delta\boldsymbol\sigma_{\rm visc}
=
2\eta
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
+
2
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}.
$$

在当前速度场 $\boldsymbol u$ 固定时，$\eta$、$\partial\eta/\partial\varepsilon_{\rm II}$ 和 $\dot{\boldsymbol\varepsilon}_{\boldsymbol u}$ 都是已知量，唯一的方向变量是

$$
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}.
$$

右端第一项对 $\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$ 是线性的。第二项中

$$
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
$$

是关于 $\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$ 的线性函数，再乘上固定的 $\dot{\boldsymbol\varepsilon}_{\boldsymbol u}$ 后仍然是线性映射。因此整体定义了一个线性映射：

$$
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\longmapsto
\delta\boldsymbol\sigma_{\rm visc}.
$$

它的输入是二阶张量，输出也是二阶张量，所以可以用四阶张量 $\mathbb C(\boldsymbol u)$ 表示：

$$
\delta\boldsymbol\sigma_{\rm visc}
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}.
$$

具体地，

$$
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
=
2\eta
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
+
2
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right)
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}.
$$

分量形式为

$$
C_{ijkl}(\boldsymbol u)
=
2\eta\,\delta_{ik}\delta_{jl}
+
2
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{ij}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{kl}.
$$

于是

$$
\left(
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right)_{ij}
=
C_{ijkl}(\boldsymbol u)
\left(\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}\right)_{kl}.
$$

把 $C_{ijkl}(\boldsymbol u)$ 代入，得到

$$
\begin{aligned}
\left(
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
\right)_{ij}
&=
2\eta\,
\delta_{ik}\delta_{jl}
\left(\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}\right)_{kl}
\\
&\quad+
2
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{ij}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{kl}
\left(\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}\right)_{kl}.
\end{aligned}
$$

第一项中，

$$
\delta_{ik}\delta_{jl}
\left(\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}\right)_{kl}
=
\left(\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}\right)_{ij}.
$$

第二项中，

$$
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{kl}
\left(\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}\right)_{kl}
=
\dot{\boldsymbol\varepsilon}_{\boldsymbol u}
:
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}.
$$

所以这个四阶张量的双点积正好还原黏性应力的一阶变化。

完整线性化应力为

$$
\delta\boldsymbol\sigma(\bm u,p)[\tilde{\boldsymbol u},\tilde p]
=
\mathbb C(\boldsymbol u):
\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}
-
\tilde p\boldsymbol I.
$$

这里 $\mathbb C(\boldsymbol u)$ 不是额外的物理模型，而是非线性黏性应力对当前速度场的一阶导数。若只保留 $2\eta\dot{\boldsymbol\varepsilon}_{\tilde{\boldsymbol u}}$，就是冻结黏度的 Picard 线性化；包含黏度导数项的完整表达式才是一致切线。

## 四阶张量的转置

四阶张量 $\mathbb C$ 的转置由二阶张量内积定义：

$$
\left(
\mathbb C:\boldsymbol A
\right):\boldsymbol B
=
\boldsymbol A:
\left(
\mathbb C^T:\boldsymbol B
\right)
\qquad
\forall \boldsymbol A,\boldsymbol B.
$$

如果

$$
\left(\mathbb C:\boldsymbol A\right)_{ij}
=
C_{ijkl}A_{kl},
$$

那么转置张量满足

$$
\left(\mathbb C^T:\boldsymbol B\right)_{kl}
=
C^T_{klij}B_{ij}.
$$

为了使

$$
\left(\mathbb C:\boldsymbol A\right):\boldsymbol B
=
\boldsymbol A:\left(\mathbb C^T:\boldsymbol B\right)
$$

对任意 $\boldsymbol A,\boldsymbol B$ 成立，需要

$$
C^T_{klij}=C_{ijkl}.
$$

等价地，若仍用 $i,j,k,l$ 作为自由指标，则

$$
C^T_{ijkl}=C_{klij}.
$$

也就是说，四阶张量转置是交换输出指标对和输入指标对：

$$
(i,j)\longleftrightarrow(k,l).
$$

这和矩阵转置完全类比。矩阵把向量映射成向量，转置交换输出指标和输入指标；这里 $\mathbb C$ 把二阶张量映射成二阶张量，所以转置交换输出指标对和输入指标对。

## 本文黏性切线的对称性

对于本文中的黏性切线，

$$
C_{ijkl}(\boldsymbol u)
=
2\eta\,\delta_{ik}\delta_{jl}
+
2
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{ij}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{kl}.
$$

交换 $(i,j)$ 和 $(k,l)$ 后，

$$
C_{klij}(\boldsymbol u)
=
2\eta\,\delta_{ki}\delta_{lj}
+
2
\frac{\partial\eta}{\partial\varepsilon_{\rm II}}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{kl}
\left(\dot{\boldsymbol\varepsilon}_{\boldsymbol u}\right)_{ij}.
$$

由于

$$
\delta_{ki}\delta_{lj}
=
\delta_{ik}\delta_{jl},
$$

且标量乘法可交换，所以

$$
C_{klij}(\boldsymbol u)
=
C_{ijkl}(\boldsymbol u).
$$

因此这个 Glen 型黏性切线满足

$$
\mathbb C(\boldsymbol u)^T
=
\mathbb C(\boldsymbol u).
$$

所以在伴随问题中，体内的伴随线性化应力与增量正问题中的线性化应力形式相同，只是方向变量从 $(\tilde{\boldsymbol u},\tilde p)$ 换成了 $(\boldsymbol v,r)$：

$$
\delta\boldsymbol\sigma^*(\bm u,p)[\boldsymbol v,r]
=
\mathbb C(\boldsymbol u)^T:
\dot{\boldsymbol\varepsilon}_{\boldsymbol v}
-
r\boldsymbol I.
$$
