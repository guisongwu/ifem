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
NonlinearStokesSlab          % 冰层正问题
NonlinearStokesMMS           % 制造解收敛测试
NonlinearStokesAdjointInversion
NonlinearStokesRegression    % 完整回归测试
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

## 正问题示例

### `NonlinearStokesSlab.m`

求解倾斜周期冰层上的非线性全 Stokes 正问题，并绘制：

- 压力；
- 水平速度；
- 竖直速度。

顶部为零牵引，底部满足不可穿透和 Weertman 滑移条件。

### `NonlinearStokesSlabContinuation.m`

冰层正问题的正则化延拓版本。依次减小
\(\varepsilon_{\mathrm{reg}}\)，并将前一级解作为下一级初值，用于提高
强非线性问题的求解稳定性。

## 制造解验证

### `NonlinearStokesMMSData.m`

定义制造解及其对应数据，包括：

- 精确速度；
- 精确压力；
- 体力；
- 顶部牵引。

### `NonlinearStokesMMS.m`

在多层网格上运行制造解测试，报告速度和压力的 \(L^2\) 误差及收敛阶。

对光滑解，P2--P1 元通常应表现为：

$$
\|u-u_h\|_{L^2}=O(h^3),
\qquad
\|p-p_h\|_{L^2}=O(h^2).
$$

### `NonlinearStokesMMSContinuation.m`

在固定网格上逐级减小正则化参数，检查 MMS 解、黏度范围和误差对
\(\varepsilon_{\mathrm{reg}}\) 的敏感性。

## 反问题

### `NonlinearStokesBetaInversion.m`

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

### `NonlinearStokesAdjointInversion.m`

伴随反演实现，采用：

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

## 诊断与测试

### `NonlinearStokesDiagnosis.m`

在相同几何、参数初值和观测条件下比较 Glen 指数 \(n=1\) 与 \(n=3\)：

- 参数 Jacobian 奇异值；
- 条件数；
- 参数误差的可观测模式；
- 最终反演误差。

### `diagnose_nonlinear_stokes_inversion_result.mat`

`NonlinearStokesDiagnosis.m` 的已保存结果，包括真实参数、初值和两组
诊断数据。

### `NonlinearStokesRegression.m`

完整 MATLAB 回归测试，覆盖：

- 三层冰层网格上的散度和完整残差；
- 压力约束的 `auto`、`none` 和 `mean-zero` 模式；
- 无牵引边界下的压力零空间处理；
- MMS 速度和压力收敛阶；
- 直接求解与正则化延拓的一致性；
- 伴随导数检查；
- 两种反演方法的稳定性和目标函数下降。

运行：

```matlab
NonlinearStokesRegression
```

全部通过时输出：

```text
All nonlinear Stokes regression checks passed.
```

该测试允许运行数分钟。

## 技术文档

### `nonlinear-stokes-tangent-matrix-explained.md`

解释残差、黏性残差、一致切线矩阵，以及 Picard 矩阵与一致切线矩阵的
区别。

### `pressure-mean-and-traction-boundary.md`

解释压力零均值规范与牵引边界之间的数学关系，以及不相容约束为什么会
导致常数散度。

### `nonlinear-stokes-residual-damping-gradient-check.md`

解释完整非线性残差、拉格朗日乘子同步阻尼和伴随梯度有限差分检查。

### 设计记录

以下文件记录开发过程中的设计决策：

- `nonlinear_stokes_continuation_design.md`
- `nonlinear-stokes-pressure-constraint-design.md`
- `nonlinear-stokes-reliability-improvements-design.md`
- `nonlinear-stokes-directory-migration-design.md`
- `matlab-file-renaming-design.md`

## 推荐验证流程

修改核心数值程序后，至少运行：

```matlab
NonlinearStokesMMS
NonlinearStokesSlab
NonlinearStokesRegression
```

重点检查：

- `info.converged` 为 `true`；
- `info.relchange(end)` 小于 `option.tol`；
- `info.nonlinearResidual(end)` 小于 `option.residual_tol`；
- MMS 速度和压力误差随网格下降；
- 散度误差随网格加密下降；
- 伴随梯度检查误差保持较小。
