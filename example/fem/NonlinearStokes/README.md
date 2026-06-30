# Nonlinear Stokes

本目录包含二维非线性全 Stokes 冰流模型的有限元求解器、正问题示例、
制造解测试、底部摩擦系数反演、诊断程序和相关技术文档。

模型采用：

- P2 速度与 P1 压力 Taylor--Hood 元；
- Glen 非线性黏度；
- Weertman 底部滑移定律；
- 左右周期边界；
- 顶部给定牵引；
- 底部不可穿透条件。

## 快速开始

在 MATLAB 中运行：

```matlab
cd('/path/to/ifem')
setpath
```

建议按以下顺序了解和验证程序：

```matlab
NSConverRate    % 制造解收敛阶测试
NSEpsContinuation
NSAdjInvTikhonov
NonlinearStokesAdjInvSlabBed
NSFDInversion
```

## 核心求解器

### `NonlinearStokesP2P1.m`

非线性 Stokes 主求解器：

```matlab
[soln,eqn,info] = NonlinearStokesP2P1(...
    node,elem,bdFlag,pde,option);
```

主要功能：

- 装配 Glen 非线性黏性项；
- 装配非线性 Weertman 底部滑移项；
- 处理周期边界和底部不可穿透约束；
- 使用阻尼 Picard 迭代；
- 检查迭代增量和完整非线性残差；
- 可选装配一致切线矩阵及参数方向导数。

主要输入数据：

```matlab
pde.A       % Glen 流动率
pde.n       % Glen 指数
pde.beta    % 底部摩擦系数，可为函数
pde.m       % Weertman 指数
pde.f       % 可选体力
pde.rho
pde.gravity
pde.g_N     % 顶部牵引
```

常用选项：

```matlab
option.periodic
option.periodic_x
option.eps_reg
option.maxIt
option.tol
option.residual_tol
option.damping
option.quadorder
option.assemble_tangent
option.pressure_constraint
```

压力约束模式包括：

```matlab
'auto'       % 有牵引边界时不固定均值，否则固定零均值
'mean-zero'  % 始终固定压力零均值
'none'       % 不添加压力均值约束
```

输出结构：

- `soln`：速度和压力；
- `eqn`：离散矩阵、约束和可选切线算子；
- `info`：迭代次数、黏度范围、完整残差和压力约束诊断。

## 制造解验证

### `NSMMSData.m`

定义制造解及其对应数据，包括：

- 精确速度；
- 精确压力；
- 体力；
- 顶部牵引。

### `NSConverRate.m`

在多层网格上运行制造解测试，报告速度和压力的 \(L^2\) 误差及收敛阶。

对光滑解，P2--P1 元通常应表现为：

$$
\|u-u_h\|_{L^2}=O(h^3),
\qquad
\|p-p_h\|_{L^2}=O(h^2).
$$

### `NSEpsContinuation.m`

在固定网格上逐级减小正则化参数，检查 MMS 解、黏度范围和误差对
\(\varepsilon_{\mathrm{reg}}\) 的敏感性。

## 反问题

### `NSFDInversion.m`

使用表面速度合成观测恢复空间变化的底部摩擦系数
\(\beta(x)\)。

算法采用：

- 参数变量 \(q=\log\beta\)，保证 \(\beta>0\)；
- 完整非线性正问题；
- 有限差分参数 Jacobian；
- Gauss--Newton/Levenberg--Marquardt 更新；
- 回溯线搜索；
- 周期一阶差分正则化。

该实现直观，适合作为反演基准，但参数维数增加后计算成本较高。

### `NSAdjInvTikhonov.m`

带周期一阶差分 Tikhonov 正则化的伴随反演实现，采用一致非线性切线矩阵、
伴随梯度和矩阵自由 Gauss--Newton 更新。

### `NonlinearStokesAdjInvSlabBed.m`

倾斜平板床上的边界积分目标函数伴随反演。算法采用：

- 一致非线性切线矩阵；
- 增量状态方程；
- 伴随梯度；
- 矩阵自由 Gauss--Newton Hessian-vector product；
- PCG 参数更新；
- 回溯线搜索。

首次迭代会检查：

- 状态方向导数；
- 伴随梯度方向导数；
- Gauss--Newton 方向量。

### `NonlinearStokesAdjInvSinBed.m`

ISMIP-HOM-B-like 正弦床上的伴随反演脚本，反演变量、观测方式和优化框架与
`NonlinearStokesAdjInvSlabBed.m` 相同。

### `NSDerivativeComparison.m`

独立比较有限差分导数与伴随导数，包括目标函数梯度和 Gauss--Newton Hessian。

## ISMIP-HOM 示例

### `ISMIPHOM_B.m`

运行二维 flowline 形式的 ISMIP-HOM experiment B，并与官方 full-Stokes 曲线和
PISM 图中提取的数据比较。

### `ISMIPHOM_D.m`

运行二维 flowline 形式的 ISMIP-HOM experiment D，并与官方 full-Stokes 曲线和
PISM 图中提取的数据比较。

### `ISMIPHOM_B_L5Fields.m` 和 `ISMIPHOM_D_L5Fields.m`

绘制 \(L=5\) km 情形下的速度和压力场。

## 归档诊断

`unused/` 中保留了较早的平板正问题、平板延拓、床面 MMS、诊断和回归脚本：

- `unused/NSSlab.m`
- `unused/NSSlabContinuation.m`
- `unused/NSConverRateBed.m`
- `unused/NSDiagnosis.m`
- `unused/NSRegression.m`

## 技术文档

### `doc/implementation/tangent-matrix.md`

解释残差、黏性残差、一致切线矩阵，以及 Picard 矩阵与一致切线矩阵的
区别。

### `doc/theory/pressure-traction.md`

解释压力零均值规范与牵引边界之间的数学关系，以及不相容约束为什么会
导致常数散度。

### `doc/implementation/residual-damping-gradient.md`

解释完整非线性残差、拉格朗日乘子同步阻尼和伴随梯度有限差分检查。

### 设计记录

以下文件记录开发过程中的设计决策：

- `doc/design/continuation.md`
- `doc/design/pressure-constraint.md`
- `doc/design/reliability.md`
- `doc/design/migration.md`
- `doc/design/renaming.md`

## 推荐验证流程

修改核心数值程序后，至少运行：

```matlab
NSConverRate
NSEpsContinuation
NSDerivativeComparison
```

重点检查：

- `info.converged` 为 `true`；
- `info.relchange(end)` 小于 `option.tol`；
- `info.nonlinearResidual(end)` 小于 `option.residual_tol`；
- MMS 速度和压力误差随网格下降；
- 散度误差随网格加密下降；
- 伴随梯度检查误差保持较小。
