# 压力零均值、牵引边界与求解器约束

## 基本关系

非线性不可压缩 Stokes 方程写作

$$
-\nabla\cdot\left(2\eta(u)\dot\varepsilon(u)\right)+\nabla p=f,
\qquad
\nabla\cdot u=0.
$$

应力和牵引为

$$
\sigma(u,p)=2\eta(u)\dot\varepsilon(u)-pI,
\qquad
t=\sigma(u,p)n.
$$

在体方程中，压力只通过 $\nabla p$ 出现。因此

$$
p\mapsto p+c
$$

不改变体内动量方程。但牵引会变成

$$
\sigma(u,p+c)n
=\sigma(u,p)n-cn.
$$

也就是

$$
\boxed{
p\mapsto p+c
\quad\Longleftrightarrow\quad
g_N\mapsto g_N-cn
}
$$

压力常数和法向牵引不能独立指定。

## 什么时候可以用压力零均值

如果整个边界只给定速度，压力常数没有被物理条件确定。这时可以加一个规范，例如

$$
\int_\Omega p\,\mathrm dx=0,
$$

它只是从等价压力族中选一个代表。

如果存在非空绝对牵引边界

$$
\sigma(u,p)n=g_N,
$$

则压力常数通常已经由 $g_N$ 固定。此时不能再任意强制压力零均值，除非牵引给出的压力本来就满足零均值。

牵引自由边界

$$
\sigma(u,p)n=0
$$

也是绝对牵引条件。它不是“没有压力条件”，而是指定总牵引为零，因此同样会固定压力参考。

## 若一定要改成零均值

设自然牵引问题得到的压力平均值为

$$
\overline p=\frac1{|\Omega|}\int_\Omega p\,\mathrm dx.
$$

若定义

$$
p_0=p-\overline p,
$$

则必须同时把牵引改为

$$
g_{N,0}=g_N+\overline p\,n.
$$

因此下面两种写法等价：

$$
\sigma(u,p)n=g_N,
\qquad
\overline p\ \text{由牵引自然确定};
$$

$$
\int_\Omega p_0\,\mathrm dx=0,
\qquad
\sigma(u,p_0)n=g_N+\overline p\,n.
$$

只平移压力而保持原来的 $g_N$ 不变，会改变边值问题。

## 离散系统中的表现

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
\end{bmatrix},
\qquad
Bu=0.
$$

压力零均值可写作

$$
m^Tp=0,
\qquad
m_i=\int_\Omega\phi_i^p\,\mathrm dx.
$$

用拉格朗日乘子 $\lambda_p$ 强制后，系统变为

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

第二行不再是 $Bu=0$，而是

$$
Bu=-m\lambda_p.
$$

如果压力零均值与牵引相容，则 $\lambda_p=0$。如果不相容，则 $\lambda_p\ne0$，离散散度会出现一个常数压力模式对应的偏移，常表现为

$$
\nabla\cdot u_h\approx\text{非零常数}.
$$

这不是普通有限元逼近误差，而是额外压力规范改变了连续性方程。

## 处理原则

1. 纯速度边界：压力常数未定，需要添加压力规范，如 $\int_\Omega p\,\mathrm dx=0$。
2. 存在绝对牵引边界：压力常数通常已由牵引固定，不应再独立强制零均值。
3. 若希望使用零均值压力表示，必须同步平移牵引：

$$
p_0=p-\overline p,
\qquad
g_{N,0}=g_N+\overline p\,n.
$$

## 求解器设计

`NonlinearStokesP2P1.m` 不应无条件添加压力均值约束。修复目标是：

- 保留顶部零牵引 `pde.g_N=[]` 的物理含义；
- 有牵引边界时默认不施加压力零均值；
- 压力确有常数零空间时继续施加零均值；
- 允许调用者显式覆盖自动判断；
- 正问题、切线系统和伴随反演使用同一个约束结构。

新增选项：

```matlab
option.pressure_constraint = 'auto';
```

允许值如下。

### `auto`

默认模式。若存在 `bdFlag == 2` 的牵引边界，则不添加压力零均值；否则添加压力零均值。

即使 `pde.g_N=[]`，`bdFlag == 2` 仍表示齐次牵引

$$
\sigma n=0,
$$

因此也属于绝对牵引边界。

### `mean-zero`

始终施加

$$
\int_\Omega p_h\,\mathrm dx=0.
$$

该模式用于制造解、特殊压力规范或兼容性实验。调用者需要保证牵引数据与零均值压力相容。

### `none`

不添加压力零均值约束。

该模式适用于压力常数已由边界条件固定的情况，也可用于诊断。若压力仍有常数零空间，线性系统可能奇异。

## 自动判断

求解器在边界识别完成后计算：

```matlab
hasTractionBoundary = any(bdFlag(:) == 2);
```

然后根据选项确定是否添加压力均值约束：

```matlab
switch option.pressure_constraint
    case 'auto'
        addPressureMeanConstraint = ~hasTractionBoundary;
    case 'mean-zero'
        addPressureMeanConstraint = true;
    case 'none'
        addPressureMeanConstraint = false;
    otherwise
        error('iFEM:NSPressureConstraint',...
            'Unknown pressure_constraint value: %s.',...
            option.pressure_constraint);
end
```

判断应使用所有 `bdFlag == 2` 的边界，而不是只依赖 `topEdgeIdx`，因为类型 2 的牵引边界不一定只在几何顶部。

## 约束矩阵

周期约束和底部不可穿透约束保持不变。只有在

```matlab
addPressureMeanConstraint == true
```

时，`buildconstraints` 才追加压力质量向量：

```matlab
pressureMean = accumarray(double(elem(:)),...
    repmat(area/3,3,1),[Np,1]);
```

因此约束矩阵 `C` 的行数会随压力约束模式变化。正问题鞍点系统仍为

$$
\begin{bmatrix}
M & C^T\\
C & 0
\end{bmatrix}.
$$

一致切线系统继续复用同一个 `C`：

$$
\begin{bmatrix}
M_{\mathrm{tan}} & C^T\\
C & 0
\end{bmatrix}.
$$

这样正问题、增量状态和伴随方程具有一致的约束空间。

## 输出信息

在 `info` 中记录最终选择，便于诊断 `'auto'` 实际采用了哪种行为：

```matlab
info.pressureConstraint = option.pressure_constraint;
info.hasTractionBoundary = hasTractionBoundary;
info.pressureMeanConstrained = addPressureMeanConstraint;
```

未知选项应立即报错。`'none'` 模式下可能出现的奇异矩阵不应被自动掩盖，因为这通常说明调用者选择了与边界条件不相容的压力规范。

## 验证计划

### `ice_slab`

默认 `'auto'` 应检测到顶部牵引边界并取消压力零均值。需要检查：

- 非线性迭代收敛；
- 鞍点矩阵满秩；
- 散度积分接近零；
- $\|\nabla\cdot u_h\|_{L^2}$ 随网格加密下降；
- 直接求解与正则化延拓结果一致。

### MMS

默认 `'auto'` 同样不固定压力均值，因为 MMS 有顶部牵引边界。需要检查：

- 速度 $L^2$ 误差约三阶；
- 压力 $L^2$ 误差约二阶；
- 散度 $L^2$ 误差随网格下降；
- 计算压力与制造解的绝对压力一致。

若显式设置

```matlab
option.pressure_constraint = 'mean-zero';
```

零均值 MMS 应与 `'auto'` 在离散误差范围内一致；而 `ice_slab` 应重现原来的常数散度，用于验证显式覆盖生效。

### `none` 和无牵引边界

对存在顶部牵引的 `ice_slab`，`'none'` 应满秩，并与 `'auto'` 一致。

对没有 `bdFlag == 2` 的测试，`'auto'` 应添加压力零均值约束，`info.pressureMeanConstrained` 应为 `true`，系统应保持可解。

### 伴随反演

运行 `NSAdjInvTikhonov.m`，检查正问题收敛、状态导数检查、伴随梯度检查和 Gauss--Newton 方向检查。

## 兼容性和风险

- 默认选项为 `'auto'`，已有调用代码不需要修改。
- 自动判断只扫描一次 `bdFlag`，不使用运行时矩阵秩检测。
- `bdFlag == 2` 被解释为绝对牵引边界；若调用者用类型 2 表示其他非标准条件，应显式设置压力约束模式。
- 某些混合边界问题可能有特殊压力零空间，自动规则不能代替完整适定性分析。
- 显式 `'mean-zero'` 允许用户构造与牵引不相容的问题，这是有意保留的高级覆盖能力。

## 决策记录

1. 保留 `ice_slab` 的顶部零牵引，不通过修改 `g_N` 修复。
2. 修改通用求解器，而不是只修改示例。
3. 默认采用基于牵引边界的自动压力规范。
4. 增加 `'mean-zero'` 和 `'none'` 显式覆盖模式。
5. 自动判断使用所有 `bdFlag == 2` 边界。
6. 不使用运行时矩阵秩检测。
7. 正问题和一致切线系统复用同一个约束矩阵。
8. 在 `info` 中公开实际采用的压力约束行为。

## 结论

对当前程序，`NonlinearStokesP2P1.m` 不应无条件添加压力均值约束。对于 `ice_slab` 这类含顶部牵引边界的问题，应去掉该约束；更一般地，应根据边界条件判断压力是否仍有常数零空间，例如：

```matlab
if pressureHasConstantNullspace
    % 添加压力零均值约束
end
```

当前边界分类下也可以先用顶部牵引边界作为简化判断：

```matlab
if isempty(topEdgeIdx)
    % 没有牵引边界时，添加压力零均值约束
end
```

核心结论：压力零均值只在压力常数未被物理边界条件确定时才是自由规范；牵引自由边界 $\sigma n=0$ 已经指定了绝对牵引，通常会固定压力常数。
