# 非线性 Stokes 方程中的残差、黏性残差与切线矩阵

本文解释 `NonlinearStokesP2P1.m` 中以下概念：

- 什么是方程残差；
- 什么是黏性残差；
- 什么是切线矩阵；
- 为什么 Glen 非线性黏度需要“一致切线矩阵”；
- Picard 矩阵与一致切线矩阵有什么区别；
- 切线矩阵为什么是伴随反演所必需的。

对应代码位于：

```text
example/fem/NonlinearStokes/NonlinearStokesP2P1.m
```

## 1. 从线性方程说起

考虑线性方程

$$
K u=f,
$$

其中：

- $u$ 是待求的离散速度；
- $K$ 是刚度矩阵；
- $f$ 是载荷向量。

把方程两边移到同一侧，可以定义残差

$$
R(u)=Ku-f.
$$

如果 $u$ 是方程的精确解，那么

$$
R(u)=0.
$$

因此，“求解方程”也可以理解为“寻找一个使残差为零的 $u$”。

对于线性问题，残差对 $u$ 的导数就是

$$
\frac{\partial R}{\partial u}=K.
$$

所以在线性问题中，刚度矩阵本身就是残差的切线矩阵。

## 2. 什么是残差

残差表示把一个候选解代入方程后，方程还差多少才能成立。

对于一般的非线性离散方程，可以写成

$$
R(U)=0,
$$

其中 $U$ 可以包含所有未知量。例如在不可压缩 Stokes 问题中，

$$
U=
\begin{bmatrix}
u\\
p\\
\lambda
\end{bmatrix},
$$

这里：

- $u$ 是速度自由度；
- $p$ 是压力自由度；
- $\lambda$ 是用于周期条件、法向约束等条件的拉格朗日乘子。

残差向量也由多个部分组成：

$$
R(U)=
\begin{bmatrix}
R_{\mathrm{momentum}}(u,p,\lambda)\\
R_{\mathrm{div}}(u)\\
R_{\mathrm{constraint}}(u)
\end{bmatrix}.
$$

分别对应：

1. 动量方程是否满足；
2. 不可压缩条件 $\nabla\cdot u=0$ 是否满足；
3. 周期边界和底部法向约束是否满足。

残差不是通常意义上的“数值误差”。它首先表示离散方程的不平衡量。残差较小通常意味着离散方程解得较准确，但它不直接等于有限元解与连续真解之间的误差。

## 3. 非线性 Stokes 动量方程

代码求解的二维非线性 Stokes 动量方程为

$$
-\nabla\cdot\left(2\eta(u)\dot{\varepsilon}(u)\right)+\nabla p=f,
$$

其中

$$
\dot{\varepsilon}(u)
=\frac{1}{2}\left(\nabla u+\nabla u^T\right)
$$

是应变率张量，$\eta(u)$ 是由 Glen 定律给出的有效黏度。

代码使用

$$
\eta(u)
=\frac{1}{2}A^{-1/n}
\left(\dot{\varepsilon}_{\mathrm{II}}(u)+\varepsilon_{\mathrm{reg}}^2\right)^{
\frac{1-n}{2n}},
$$

其中

$$
\dot{\varepsilon}_{\mathrm{II}}(u)
=\frac{1}{2}\dot{\varepsilon}(u):\dot{\varepsilon}(u).
$$

符号“$:$”表示张量双点积。在二维情况下，

$$
\dot{\varepsilon}:\dot{\varepsilon}
=\dot{\varepsilon}_{xx}^2+\dot{\varepsilon}_{zz}^2
+2\dot{\varepsilon}_{xz}^2.
$$

当 $n\ne 1$ 时，黏度依赖于应变率，而应变率又依赖于速度。因此动量方程是关于速度 $u$ 的非线性方程。

正则化参数 $\varepsilon_{\mathrm{reg}}>0$ 用于避免零应变率附近黏度奇异或导数无界。

## 4. 什么是黏性残差

有限元方法不直接在每个空间点强制满足微分方程，而是要求方程对所有测试函数 $v$ 满足弱形式。

动量方程中的黏性项为

$$
-\nabla\cdot\left(2\eta(u)\dot{\varepsilon}(u)\right).
$$

经过分部积分后，其弱形式主要部分为

$$
r_{\mathrm{visc}}(u;v)
=\int_\Omega
2\eta(u)\,
\dot{\varepsilon}(u):\dot{\varepsilon}(v)\,\mathrm{d}x.
$$

这就是黏性残差在测试函数 $v$ 方向上的作用。

将所有有限元基函数依次作为测试函数，就得到一个离散向量

$$
R_{\mathrm{visc}}(u).
$$

它表示当前速度 $u$ 产生的内部黏性力。完整的动量残差还应包含压力、体力、边界摩擦和边界外力等部分。例如可概括为

$$
R_{\mathrm{momentum}}
=R_{\mathrm{visc}}(u)
+R_{\mathrm{bed}}(u,\beta)
+B^Tp
+C^T\lambda
-F.
$$

所以：

- `黏性残差`只是完整残差中由冰体内部黏性应力产生的部分；
- `完整残差`还包括压力、不可压缩条件、底部滑移、载荷和约束。

## 5. 什么是切线矩阵

对非线性残差

$$
R(U)=0
$$

在当前状态 $U$ 附近施加一个小扰动 $\delta U$。一阶 Taylor 展开为

$$
R(U+\delta U)
\approx
R(U)+R_U(U)\,\delta U,
$$

其中

$$
R_U(U)=\frac{\partial R}{\partial U}
$$

称为：

- 切线矩阵；
- Jacobian 矩阵；
- 一致线性化矩阵；
- Newton 矩阵。

这些名称在本文语境中指的是同一个对象。

它描述当前状态附近“未知量改变一点，残差会怎样改变”：

$$
\delta R\approx R_U\,\delta U.
$$

“切线”一词来自一维函数。对标量函数 $r(u)$，

$$
r(u+\delta u)\approx r(u)+r'(u)\delta u.
$$

$r'(u)$ 是曲线在当前点的切线斜率。多维方程中的 $R_U$ 就是这个斜率的矩阵推广。

## 6. 黏性残差的一致线性化

黏性弱形式为

$$
r_{\mathrm{visc}}(u;v)
=\int_\Omega
2\eta(u)\,
\dot{\varepsilon}(u):\dot{\varepsilon}(v)\,\mathrm{d}x.
$$

考虑速度方向 $\delta u$。因为 $\eta$ 和 $\dot{\varepsilon}(u)$ 都依赖于 $u$，使用乘积求导可得

$$
\begin{aligned}
D r_{\mathrm{visc}}(u)[\delta u;v]
= {}&
\int_\Omega
2\eta(u)\,
\dot{\varepsilon}(\delta u):\dot{\varepsilon}(v)\,\mathrm{d}x\\
&+
\int_\Omega
2D\eta(u)[\delta u]\,
\dot{\varepsilon}(u):\dot{\varepsilon}(v)\,\mathrm{d}x.
\end{aligned}
$$

第一项来自应变率 $\dot{\varepsilon}(u)$ 的变化，第二项来自黏度 $\eta(u)$ 的变化。

### 6.1 冻结黏度项

第一项是

$$
\int_\Omega
2\eta(u)\,
\dot{\varepsilon}(\delta u):\dot{\varepsilon}(v)\,\mathrm{d}x.
$$

如果把当前黏度 $\eta(u)$ 看作一个已知系数，只装配这一项，就得到代码中 Picard 迭代使用的黏性矩阵。

### 6.2 黏度导数项

令

$$
s=\dot{\varepsilon}_{\mathrm{II}}(u)+\varepsilon_{\mathrm{reg}}^2,
\qquad
a=\frac{1-n}{2n}.
$$

则

$$
\eta(u)=\frac12 A^{-1/n}s^a.
$$

对速度方向 $\delta u$ 求导：

$$
D\eta(u)[\delta u]
=\eta(u)\frac{a}{s}
\left(\dot{\varepsilon}(u):\dot{\varepsilon}(\delta u)\right).
$$

代回黏性残差的导数后，得到

$$
\begin{aligned}
D r_{\mathrm{visc}}(u)[\delta u;v]
=\int_\Omega 2\eta(u)\Bigg[
&\dot{\varepsilon}(\delta u):\dot{\varepsilon}(v)\\
&+\frac{a}{s}
\big(\dot{\varepsilon}(u):\dot{\varepsilon}(\delta u)\big)
\big(\dot{\varepsilon}(u):\dot{\varepsilon}(v)\big)
\Bigg]\mathrm{d}x.
\end{aligned}
$$

`assembleviscoustangent` 装配的就是这个完整表达式。

因为它是从原非线性黏性残差直接求导得到的，所以称为“一致切线矩阵”。

## 7. Picard 矩阵与一致切线矩阵

两者的差别可以概括为：

| 矩阵 | 是否考虑 $\dot{\varepsilon}(u)$ 的变化 | 是否考虑 $\eta(u)$ 的变化 | 主要用途 |
|---|---:|---:|---|
| Picard 矩阵 | 是 | 否 | 稳健地求非线性正问题 |
| 一致切线矩阵 | 是 | 是 | Newton 法、灵敏度、伴随和反演 |

Picard 方法在第 $k$ 次迭代中先计算

$$
\eta^k=\eta(u^k),
$$

然后冻结这个黏度，求解一个线性 Stokes 问题：

$$
-\nabla\cdot\left(2\eta^k\dot{\varepsilon}(u^{k+1})\right)
+\nabla p^{k+1}=f.
$$

它没有声称冻结黏度矩阵就是原非线性方程的完整导数。它只是构造一个便于迭代的线性问题。

一致切线矩阵则必须回答：

> 如果当前速度发生一个微小变化，包含黏度变化在内的完整残差会怎样变化？

因此，一致切线必须包含 $D\eta(u)[\delta u]$。

## 8. 底部滑移残差及其切线矩阵

底部 Weertman 滑移关系可写为

$$
t_b(u_t,\beta)
=\beta
\left(u_t^2+\varepsilon_{\mathrm{reg}}^2\right)^{(m-1)/2}
u_t,
$$

其中：

- $u_t$ 是底部切向速度；
- $\beta$ 是待反演的底部摩擦系数；
- $m$ 是滑移定律指数。

相应的底部弱形式残差为

$$
r_{\mathrm{bed}}(u;v)
=\int_{\Gamma_b}t_b(u_t,\beta)v_t\,\mathrm{d}s.
$$

令

$$
s_b=u_t^2+\varepsilon_{\mathrm{reg}}^2.
$$

底部牵引对切向速度的导数为

$$
\frac{\partial t_b}{\partial u_t}
=\beta\left[
s_b^{(m-1)/2}
+(m-1)u_t^2s_b^{(m-3)/2}
\right].
$$

因此底部残差在速度方向 $\delta u$ 上的导数为

$$
D r_{\mathrm{bed}}(u)[\delta u;v]
=\int_{\Gamma_b}
\frac{\partial t_b}{\partial u_t}
\delta u_t\,v_t\,\mathrm{d}s.
$$

`assemblebedtangent` 装配的就是这一项。

当 $m=1$ 时，

$$
t_b=\beta u_t,
\qquad
\frac{\partial t_b}{\partial u_t}=\beta.
$$

此时底部滑移项关于速度是线性的，普通底部矩阵和底部切线矩阵相同。

当 $m\ne1$ 时，底部牵引随速度非线性变化，切线矩阵必须包含额外的导数项。

## 9. 完整 Stokes 切线系统

把速度、压力和约束放在一起，代码最终构造的切线系统具有如下分块形式：

$$
R_U=
\begin{bmatrix}
K_{\mathrm{visc}}^{\mathrm{tan}}
+K_{\mathrm{bed}}^{\mathrm{tan}} & B^T & C^T\\
B & 0 & 0\\
C & 0 & 0
\end{bmatrix}.
$$

其中：

- $K_{\mathrm{visc}}^{\mathrm{tan}}$ 是 Glen 黏性残差的一致切线；
- $K_{\mathrm{bed}}^{\mathrm{tan}}$ 是底部滑移残差的一致切线；
- $B$ 是速度散度矩阵；
- $C$ 表示周期条件、底部不可穿透条件和压力规范化等约束。

代码将这个矩阵保存在

```matlab
eqn.tangent
```

中。只有设置

```matlab
option.assemble_tangent = true;
```

时才会额外装配它。

普通正问题默认不装配切线矩阵，因为 Picard 迭代已经可以完成正问题求解，而一致切线会增加装配时间和内存开销。

## 10. 为什么反演需要一致切线矩阵

反问题中的状态方程可写为

$$
R(U,\beta)=0,
$$

其中 $\beta$ 是底部摩擦参数。

现在给参数一个小扰动 $\delta\beta$。状态也会产生相应变化 $\delta U$。对状态方程求导：

$$
R_U\,\delta U+R_\beta\,\delta\beta=0.
$$

所以增量状态满足

$$
R_U\,\delta U=-R_\beta\,\delta\beta.
$$

这里：

- `eqn.tangent` 对应 $R_U$；
- `eqn.applyBetaDerivative` 计算 $R_\beta\delta\beta$；
- $\delta U$ 是参数扰动导致的速度、压力和约束变量变化。

如果用 Picard 矩阵代替 $R_U$，就会漏掉黏度随速度变化的导数；当底部滑移是非线性的，还会漏掉滑移牵引随速度变化的导数。由此得到的 $\delta U$ 不再是原非线性方程的正确一阶灵敏度。

## 11. 切线矩阵与伴随方程

设目标函数为

$$
J(U,\beta),
$$

例如衡量计算得到的表面速度与观测速度之间的差异。

目标函数对参数的总导数包含状态变化：

$$
\frac{\mathrm{d}J}{\mathrm{d}\beta}
=J_\beta+J_U\frac{\mathrm{d}U}{\mathrm{d}\beta}.
$$

直接计算 $\mathrm{d}U/\mathrm{d}\beta$ 需要为许多参数方向求解增量状态方程。伴随法引入伴随变量 $\psi$，满足

$$
R_U^T\psi=-J_U^T.
$$

于是梯度可以写为

$$
\frac{\mathrm{d}J}{\mathrm{d}\beta}
=J_\beta+\psi^T R_\beta.
$$

`NSAdjointInversion.m` 中的核心代码为

```matlab
adjoint = eqn.tangent'\(-observationGradient);
gradient = G'*adjoint + regularizationGradient;
```

其中：

- `eqn.tangent'` 是一致切线矩阵的转置；
- `adjoint` 是伴随变量；
- `G` 是参数导数矩阵，即离散形式的 $R_\beta$；
- `gradient` 是目标函数对参数的梯度。

因此伴随法并不是只需要“某个能解正问题的矩阵”，而是需要原非线性残差的真实 Jacobian。

## 12. `applyBetaDerivative` 的含义

底部牵引为

$$
t_b(u_t,\beta)
=\beta s_b^{(m-1)/2}u_t.
$$

保持当前状态 $u$ 不变，只让参数产生扰动 $\delta\beta$，则

$$
D_\beta t_b[\delta\beta]
=\delta\beta\,s_b^{(m-1)/2}u_t.
$$

把这一方向导数投影到有限元测试函数上，就得到

$$
R_\beta\,\delta\beta.
$$

代码中的

```matlab
eqn.applyBetaDerivative(betaDirection)
```

计算的就是这个向量。它不是一个新的正问题解，而是参数变化对当前方程残差产生的一阶影响。

## 13. 一个标量类比

考虑非线性方程

$$
R(u,\beta)=\beta u^3-f=0.
$$

它对状态的切线为

$$
R_u=3\beta u^2,
$$

对参数方向 $\delta\beta$ 的导数为

$$
R_\beta\delta\beta=u^3\delta\beta.
$$

参数改变后，状态变化满足

$$
3\beta u^2\delta u=-u^3\delta\beta.
$$

所以

$$
\delta u=-\frac{u}{3\beta}\delta\beta.
$$

如果为了迭代方便，把 $u^2$ 冻结并构造另一个近似矩阵，它可能仍能用于求解原方程；但只有真实导数 $3\beta u^2$ 才能给出正确的一阶灵敏度。

非线性 Stokes 中的一致切线矩阵就是这个标量导数在有限元、多变量和约束系统中的推广。

## 14. 如何验证切线矩阵是否正确

给定一个状态方向 $\delta U$ 和小参数 $h$，可以比较：

$$
\frac{R(U+h\delta U)-R(U-h\delta U)}{2h}
$$

与

$$
R_U(U)\delta U.
$$

若切线实现正确，并且 $h$ 处于合适范围，两者的相对误差应当较小：

$$
\frac{
\left\|
\dfrac{R(U+h\delta U)-R(U-h\delta U)}{2h}
-R_U\delta U
\right\|
}{
\left\|
\dfrac{R(U+h\delta U)-R(U-h\delta U)}{2h}
\right\|
}
\ll1.
$$

反演代码还通过比较目标函数有限差分与伴随梯度方向导数来验证整条导数链：

$$
\frac{J(\beta+h\delta\beta)-J(\beta-h\delta\beta)}{2h}
\approx
\nabla_\beta J\cdot\delta\beta.
$$

这种检查可以同时发现以下问题：

- 黏性切线漏项或系数错误；
- 底部切线错误；
- 参数导数 $R_\beta$ 错误；
- 伴随方程符号错误；
- 观测算子或参数映射实现错误。

## 15. 常见误解

### 15.1 残差不等于解误差

残差衡量离散方程是否满足；解误差衡量离散解与连续真解之间的差别。二者有关，但不是同一个量。

### 15.2 切线矩阵不是几何网格的切线

这里的“切线”指非线性残差在当前解附近的一阶导数，与边界曲线或网格边的几何切向量不是同一概念。

代码中底部滑移项同时使用了底部几何切向量和非线性切线矩阵，两者不要混淆：

- 几何切向量用于提取底部切向速度 $u_t$；
- 切线矩阵描述残差对所有状态自由度的导数。

### 15.3 Picard 矩阵不是错误矩阵

Picard 矩阵适合构造稳健的非线性正问题迭代。问题只在于，它通常不是完整非线性残差的 Jacobian，因此不能不加分析地用于精确灵敏度和伴随梯度。

### 15.4 正问题收敛不代表导数一定正确

即使 Picard 迭代得到非常小的相对变化量，也只能说明正问题迭代收敛。切线矩阵、参数导数和伴随梯度仍需要单独进行有限差分检查。

## 16. 总结

可以用四句话概括：

1. 残差 $R(U)$ 表示当前候选解代入离散方程后还剩下多少不平衡。
2. 黏性残差是完整 Stokes 残差中由内部黏性应力产生的部分。
3. 切线矩阵 $R_U$ 是残差对状态变量的 Jacobian，描述状态微小变化引起的残差一阶变化。
4. Glen 黏度和非线性底部滑移都依赖速度，因此一致切线必须同时包含应变率变化、黏度变化和滑移牵引变化产生的导数项。

正问题的 Picard 迭代可以使用冻结系数的近似矩阵；非线性灵敏度、伴随梯度和 Gauss--Newton 反演则需要与原残差一致的切线矩阵。
