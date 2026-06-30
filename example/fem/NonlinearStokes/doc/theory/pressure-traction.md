# Stokes 方程中压力零均值与牵引边界的关系

本文说明不可压缩 Stokes 方程中以下问题：

- 为什么压力有时只能确定到一个常数；
- 压力零均值条件什么时候可以使用；
- 顶部牵引条件为什么会固定压力常数；
- 为什么不能同时任意指定压力零均值和绝对牵引；
- `unused/NSSlab.m` 中为什么出现非零常数散度；
- MMS 测试为什么没有出现同样的问题。

相关代码：

```text
example/fem/NonlinearStokes/NonlinearStokesP2P1.m
example/fem/NonlinearStokes/unused/NSSlab.m
example/fem/NonlinearStokes/NSConverRate.m
example/fem/NonlinearStokes/NSMMSData.m
```

## 1. 不可压缩 Stokes 方程

考虑非线性不可压缩 Stokes 方程：

$$
-\nabla\cdot\left(2\eta(u)\dot{\varepsilon}(u)\right)
+\nabla p=f
\qquad\text{in }\Omega,
$$

$$
\nabla\cdot u=0
\qquad\text{in }\Omega.
$$

其中：

- $u$ 是速度；
- $p$ 是压力；
- $\eta(u)$ 是有效黏度；
- $\dot{\varepsilon}(u)$ 是应变率张量，

$$
\dot{\varepsilon}(u)
=\frac12\left(\nabla u+\nabla u^T\right).
$$

Cauchy 应力张量为

$$
\sigma(u,p)
=2\eta(u)\dot{\varepsilon}(u)-pI.
$$

边界上的牵引为

$$
t=\sigma(u,p)n,
$$

其中 $n$ 是外法向量。

## 2. 为什么压力有时只能确定到一个常数

在区域内部，压力只通过梯度 $\nabla p$ 出现在动量方程中。

如果把压力替换为

$$
p_c=p+c,
$$

其中 $c$ 是任意常数，那么

$$
\nabla p_c=\nabla(p+c)=\nabla p.
$$

所以压力加常数不会改变区域内部的动量方程。

如果整个边界只给定速度，例如

$$
u=g_D
\qquad\text{on }\partial\Omega,
$$

那么边界条件也不直接涉及压力。因此，如果 $(u,p)$ 是解，则 $(u,p+c)$ 也是解。

这种情况下，压力只确定到一个加法常数。为了选出唯一压力，可以附加一种压力规范，例如：

$$
\int_\Omega p\,\mathrm{d}x=0,
$$

或者固定某一点的压力：

$$
p(x_0)=0.
$$

压力零均值在这里不是额外的物理条件，而只是从一族等价压力中选择一个代表。

## 3. 压力常数会改变边界牵引

虽然压力加常数不改变区域内部的 $\nabla p$，但会改变应力和边界牵引。

令

$$
p_c=p+c.
$$

对应的应力为

$$
\begin{aligned}
\sigma(u,p_c)
&=2\eta\dot{\varepsilon}(u)-(p+c)I\\
&=\sigma(u,p)-cI.
\end{aligned}
$$

因此边界牵引变成

$$
\begin{aligned}
t_c
&=\sigma(u,p_c)n\\
&=\sigma(u,p)n-cn\\
&=t-cn.
\end{aligned}
$$

所以有：

$$
\boxed{
p\mapsto p+c
\quad\Longrightarrow\quad
t\mapsto t-cn
}
$$

这说明压力常数与边界法向牵引直接相关。

## 4. 顶部牵引条件会固定压力常数

如果顶部边界 $\Gamma_t$ 给定绝对牵引

$$
\sigma(u,p)n=g_N
\qquad\text{on }\Gamma_t,
$$

那么把压力改为 $p+c$ 后，左侧变为

$$
\sigma(u,p+c)n=g_N-cn.
$$

只要给定的 $g_N$ 保持不变，通常必须有

$$
c=0.
$$

因此，非空牵引边界通常已经确定了压力常数。

这时压力不再具有任意加常数的自由度，也就不能再随意附加

$$
\int_\Omega p\,\mathrm{d}x=0.
$$

除非给定牵引数据恰好与零均值压力规范相容。

## 5. 牵引自由边界并不表示压力任意

在冰层顶部常使用牵引自由条件：

$$
\sigma(u,p)n=0.
$$

它表示顶部总牵引为零，并不表示没有压力条件。

将应力展开：

$$
2\eta\dot{\varepsilon}(u)n-pn=0.
$$

因此顶部压力与黏性法向应力必须相互平衡：

$$
pn=2\eta\dot{\varepsilon}(u)n.
$$

在近似静水状态下，顶部黏性法向应力较小，于是顶部压力接近零。

所以牵引自由边界实际上给出了压力绝对参考，而不是留下任意压力常数。

## 6. 静水压力例子

考虑高度为 $H=1$ 的冰层，使用竖直坐标 $q\in[0,1]$，顶部为 $q=1$，底部为 $q=0$。

在重力作用下，静水压力近似为

$$
p(q)=1-q.
$$

顶部压力为

$$
p(1)=0,
$$

符合顶部牵引自由条件。

但是压力平均值为

$$
\overline p
=\int_0^1(1-q)\,\mathrm{d}q
=\frac12.
$$

因此一个物理上合理的牵引自由解可以同时满足

$$
p_{\mathrm{top}}=0,
\qquad
\overline p=\frac12,
$$

但不能同时满足

$$
p_{\mathrm{top}}=0,
\qquad
\overline p=0.
$$

## 7. 如果一定要把压力改成零均值

设自然压力的平均值为

$$
\overline p
=\frac{1}{|\Omega|}
\int_\Omega p\,\mathrm{d}x.
$$

定义零均值压力

$$
p_0=p-\overline p.
$$

则

$$
\int_\Omega p_0\,\mathrm{d}x=0.
$$

对应应力为

$$
\begin{aligned}
\sigma(u,p_0)
&=2\eta\dot{\varepsilon}(u)-(p-\overline p)I\\
&=\sigma(u,p)+\overline p I.
\end{aligned}
$$

所以边界牵引必须同时改为

$$
\sigma(u,p_0)n
=g_N+\overline p\,n.
$$

因此，以下两组条件描述相同的速度场：

### 自然压力规范

$$
\sigma(u,p)n=g_N,
\qquad
\overline p\ \text{由方程和牵引自然确定}.
$$

### 零均值压力规范

$$
\int_\Omega p_0\,\mathrm{d}x=0,
$$

$$
\sigma(u,p_0)n
=g_N+\overline p\,n.
$$

不能只做

$$
p\mapsto p-\overline p
$$

而保持原来的 $g_N$ 不变。那会改变边值问题。

## 8. 弱形式中的压力和牵引

动量方程乘测试函数 $v$ 并分部积分，得到

$$
\int_\Omega
2\eta\dot{\varepsilon}(u):\dot{\varepsilon}(v)\,\mathrm{d}x
-\int_\Omega p\,\nabla\cdot v\,\mathrm{d}x
=\int_\Omega f\cdot v\,\mathrm{d}x
+\int_{\Gamma_t}g_N\cdot v\,\mathrm{d}s.
$$

将压力改为 $p+c$ 后，压力项改变为

$$
-\int_\Omega (p+c)\nabla\cdot v\,\mathrm{d}x
=-\int_\Omega p\nabla\cdot v\,\mathrm{d}x
-c\int_\Omega\nabla\cdot v\,\mathrm{d}x.
$$

根据散度定理，

$$
\int_\Omega\nabla\cdot v\,\mathrm{d}x
=\int_{\partial\Omega}v\cdot n\,\mathrm{d}s.
$$

因此额外项为

$$
-c\int_{\partial\Omega}v\cdot n\,\mathrm{d}s.
$$

它正好等价于把边界牵引增加

$$
-cn.
$$

这从弱形式再次说明：当测试函数在牵引边界上可以有非零法向分量时，压力常数不能与牵引数据相互独立地指定。

## 9. 有限元离散中的压力零均值约束

离散 Stokes 系统可写成

$$
\begin{bmatrix}
K & B^T\\
B & 0
\end{bmatrix}
\begin{bmatrix}
u\\
p
\end{bmatrix}
=
\begin{bmatrix}
F\\
0
\end{bmatrix}.
$$

其中：

- $K$ 是黏性和底部滑移矩阵；
- $B$ 是离散散度矩阵；
- $B^T$ 是离散压力梯度矩阵。

不可压缩方程是

$$
Bu=0.
$$

压力零均值可以写成

$$
m^Tp=0,
$$

其中 $m$ 是压力质量向量：

$$
m_i=\int_\Omega\phi_i^p\,\mathrm{d}x.
$$

用拉格朗日乘子 $\lambda_p$ 强制这个条件后，系统变成

$$
\begin{bmatrix}
K & B^T & 0\\
B & 0 & m\\
0 & m^T & 0
\end{bmatrix}
\begin{bmatrix}
u\\
p\\
\lambda_p
\end{bmatrix}
=
\begin{bmatrix}
F\\
0\\
0
\end{bmatrix}.
$$

其中连续性方程不再是

$$
Bu=0,
$$

而是

$$
Bu+m\lambda_p=0.
$$

如果压力零均值与牵引条件兼容，则

$$
\lambda_p=0
$$

并仍有

$$
Bu=0.
$$

如果二者不兼容，则

$$
\lambda_p\ne0,
$$

于是

$$
Bu=-m\lambda_p.
$$

因为 $m$ 代表常数压力测试模式，这通常表现为

$$
\nabla\cdot u_h=\text{非零常数}.
$$

这不是普通的网格离散误差，而是增广约束改变了连续性方程。

## 10. `ice_slab` 中出现了什么

`unused/NSSlab.m` 使用：

$$
g_N=0
$$

作为顶部牵引自由条件，同时 `NonlinearStokesP2P1.m` 无条件加入了

$$
\int_\Omega p_h\,\mathrm{d}x=0.
$$

重力产生的自然压力近似为静水压力，其平均值约为

$$
\overline p\approx0.5.
$$

因此顶部零牵引和压力零均值不兼容。

MATLAB 数值检查得到：

$$
\|\nabla\cdot u_h\|_{L^2}
\approx6.228\times10^{-2}.
$$

散度在区域内几乎为常数：

```text
min(div u_h)  = -6.2280158e-2
max(div u_h)  = -6.2279321e-2
mean(div u_h) = -6.2279724e-2
```

并且该数值在三层网格上不下降：

```text
h = 1/8  : 6.22797243e-2
h = 1/16 : 6.22797255e-2
h = 1/32 : 6.22797255e-2
```

这说明它来自边界条件和压力规范的不兼容，而不是有限元逼近误差。

## 11. 顶部牵引平移的数值验证

保留压力零均值约束，把顶部牵引改为

$$
g_N=c\,n.
$$

数值结果为：

| $c$ | $\|\nabla\cdot u_h\|_{L^2}$ |
|---:|---:|
| 0 | $6.2280\times10^{-2}$ |
| 0.25 | $8.1595\times10^{-3}$ |
| 0.50 | $1.6864\times10^{-5}$ |

当

$$
c\approx0.5
$$

时散度几乎消失。

这与静水压力的自然平均值

$$
\overline p\approx0.5
$$

一致，验证了压力平移和牵引平移之间的数学关系。

但这不应作为 `ice_slab` 的最终修复，因为设置

$$
g_N=0.5n
$$

改变了原来“顶部牵引为零”的物理边界条件。它只是一个用于验证问题来源的对照实验。

## 12. MMS 为什么没有出现问题

MMS 使用的精确压力为

$$
p(x,q)
=\cos(2\pi x/L)(1+q).
$$

由于余弦在完整周期上的平均值为零，

$$
\int_0^L\cos(2\pi x/L)\,\mathrm{d}x=0,
$$

所以

$$
\int_\Omega p\,\mathrm{d}x=0.
$$

MMS 的顶部牵引又是根据这个精确压力和精确速度计算得到的：

$$
g_N=\sigma(u,p)n.
$$

因此，MMS 的压力零均值与顶部牵引数据是相容的。

数值测试得到：

| $h$ | $\|\nabla\cdot u_h\|_{L^2}$ | 收敛阶 |
|---:|---:|---:|
| 1/4 | $4.7982\times10^{-2}$ | — |
| 1/8 | $1.4490\times10^{-2}$ | 1.73 |
| 1/16 | $4.0019\times10^{-3}$ | 1.86 |

同时散度积分接近机器零：

$$
\int_\Omega\nabla\cdot u_h\,\mathrm{d}x
\approx10^{-15}.
$$

所以 MMS 中的散度是正常有限元误差，并随网格加密趋于零。

## 13. 正确的处理原则

### 情形一：整个边界给定速度

压力只确定到常数，需要选择压力规范：

$$
\int_\Omega p\,\mathrm{d}x=0
$$

或

$$
p(x_0)=0.
$$

### 情形二：存在绝对牵引边界

如果

$$
\sigma(u,p)n=g_N
$$

在非空边界上给定，压力常数通常已经由牵引固定。

此时不应再独立强制压力零均值。

### 情形三：希望使用零均值压力表示

可以把压力平移为零均值，但必须同时平移牵引：

$$
p_0=p-\overline p,
$$

$$
g_{N,0}=g_N+\overline p\,n.
$$

这只是同一个物理解的不同压力规范。

## 14. 对当前程序的建议

当前 `NonlinearStokesP2P1.m` 无条件添加压力均值约束：

```matlab
pressureMean = accumarray(double(elem(:)),...
    repmat(area/3,3,1),[Np,1]);
...
C = sparse(double(I),double(J),S,row,2*Nu+Np);
```

对于存在顶部牵引边界的 `ice_slab`，应去掉该压力均值约束。

数值检查表明，删除它之后的鞍点矩阵仍然满秩：

```text
删除压力均值约束后的矩阵阶数：718
结构秩：718
```

这说明顶部牵引已经固定了压力常数，不需要额外规范化。

更一般的实现应根据边界条件决定是否添加压力规范，而不是无条件添加。例如：

```matlab
if pressureHasConstantNullspace
    % 添加压力零均值约束
end
```

对于当前边界分类，可以首先采用：

```matlab
if isempty(topEdgeIdx)
    % 没有牵引边界时，添加压力零均值约束
end
```

实际通用实现还应考虑其他牵引边界、混合边界以及用户是否显式要求某种压力规范。

## 15. 总结

核心关系是：

$$
\boxed{
p\mapsto p+c
\quad\Longleftrightarrow\quad
g_N\mapsto g_N-cn
}
$$

因此：

1. 压力零均值只在压力常数未被物理边界条件确定时才是自由规范；
2. 给定绝对牵引通常会固定压力常数；
3. 牵引自由边界 $\sigma n=0$ 也是绝对牵引条件；
4. 如果平移压力，就必须同步平移法向牵引；
5. `ice_slab` 同时强制顶部零牵引和压力零均值，导致离散系统产生非零常数散度；
6. 正确修复是保留物理牵引条件，并删除多余的压力零均值约束。
