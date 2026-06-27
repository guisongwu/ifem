# 非线性 Stokes 压力约束修复设计

## 1. 背景

`NonlinearStokesP2P1.m` 当前无条件施加压力零均值约束：

$$
\int_\Omega p_h\,\mathrm{d}x=0.
$$

当边界上存在给定牵引

$$
\sigma(u,p)n=g_N
$$

时，牵引条件通常已经确定压力常数。此时再强制压力零均值会过度约束压力，并通过压力均值约束的拉格朗日乘子把不可压缩方程改变为

$$
Bu+m\lambda_p=0.
$$

在 `NSSlab.m` 中，这表现为不随网格加密下降的非零常数散度。

## 2. 目标

- 修复通用求解器，而不只修补 `ice_slab` 示例。
- 保留顶部零牵引 `pde.g_N=[]` 的物理含义。
- 有牵引边界时默认不施加压力零均值约束。
- 压力确有常数零空间时继续提供零均值约束。
- 允许调用者显式覆盖自动判断。
- 保持现有调用方式兼容。
- 使正问题、切线系统、伴随反演使用相同的约束结构。

## 3. 非目标

- 不改变 Glen 黏度模型。
- 不改变 Weertman 底部滑移模型。
- 不改变周期边界或底部不可穿透约束。
- 不通过修改顶部牵引数据补偿压力平移。
- 不在每次求解时用数值秩检测决定压力规范。

## 4. 新增选项

新增求解器选项：

```matlab
option.pressure_constraint = 'auto';
```

允许值如下。

### `auto`

默认模式：

```text
存在 bdFlag == 2 的牵引边界 → 不添加压力零均值约束
不存在牵引边界             → 添加压力零均值约束
```

即使 `pde.g_N=[]`，`bdFlag==2` 仍表示齐次牵引

$$
\sigma n=0,
$$

因此仍属于绝对牵引边界。

### `mean-zero`

无论边界类型如何，始终施加

$$
\int_\Omega p_h\,\mathrm{d}x=0.
$$

该模式用于制造解、特殊压力规范或兼容性实验。调用者有责任确保给定牵引与零均值压力相容。

### `none`

不添加压力零均值约束。

该模式适用于压力常数已由边界条件固定的情况，也可用于诊断。若压力仍有常数零空间，线性系统可能奇异。

## 5. 自动判断

求解器在边界边识别完成后计算：

```matlab
hasTractionBoundary = any(bdFlag(:) == 2);
```

根据选项确定：

```matlab
switch option.pressure_constraint
    case 'auto'
        addPressureMeanConstraint = ~hasTractionBoundary;
    case 'mean-zero'
        addPressureMeanConstraint = true;
    case 'none'
        addPressureMeanConstraint = false;
    otherwise
        error(...);
end
```

判断使用所有 `bdFlag==2` 的边界，而不是只依赖变量名 `topEdgeIdx`，因为类型 2 的牵引边界不一定只位于几何顶部。

## 6. 约束矩阵

周期约束和底部不可穿透约束保持不变。

`buildconstraints` 仅在

```matlab
addPressureMeanConstraint == true
```

时追加压力质量向量：

```matlab
pressureMean = accumarray(double(elem(:)),...
    repmat(area/3,3,1),[Np,1]);
```

因此约束矩阵 `C` 的行数随压力约束模式变化。

正问题鞍点系统仍为

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

这样正问题、增量状态和伴随方程具有完全一致的约束空间。

## 7. 输出和可诊断性

在 `info` 中记录最终选择：

```matlab
info.pressureConstraint = option.pressure_constraint;
info.hasTractionBoundary = hasTractionBoundary;
info.pressureMeanConstrained = addPressureMeanConstraint;
```

这使用户可以确认 `'auto'` 实际选择了哪种行为。

## 8. 错误处理

对于未知选项值，立即报错：

```matlab
error('iFEM:NSPressureConstraint',...
    'Unknown pressure_constraint value: %s.',...);
```

不自动捕获或掩盖 `'none'` 模式下可能出现的奇异矩阵，因为这通常表示调用者选择了与边界条件不相容的压力规范。

## 9. 验证计划

### 9.1 `ice_slab`

默认 `'auto'` 应检测到顶部牵引边界并取消压力零均值。

检查：

- 非线性迭代收敛；
- 鞍点矩阵满秩；
- 散度积分接近零；
- $\|\nabla\cdot u_h\|_{L^2}$ 随网格加密下降；
- 直接求解与正则化延拓结果一致。

### 9.2 MMS

默认 `'auto'` 同样不固定压力均值，因为 MMS 有顶部牵引边界。

检查：

- 速度 $L^2$ 误差约三阶；
- 压力 $L^2$ 误差约二阶；
- 散度 $L^2$ 误差随网格下降；
- 计算压力与制造解的绝对压力一致。

### 9.3 显式 `mean-zero`

对零均值 MMS：

```matlab
option.pressure_constraint = 'mean-zero';
```

结果应与 `'auto'` 在离散误差范围内一致。

对 `ice_slab`，该模式预计重现原来的常数散度，用于证明显式覆盖确实生效。

### 9.4 显式 `none`

对存在顶部牵引的 `ice_slab`：

- 系统应满秩；
- 结果应与 `'auto'` 一致。

### 9.5 无牵引边界

构造没有 `bdFlag==2` 的测试：

- `'auto'` 应添加压力零均值约束；
- `info.pressureMeanConstrained` 应为 `true`；
- 系统应保持可解。

### 9.6 伴随反演

运行 `NSAdjointInversion.m`：

- 正问题收敛；
- 状态导数检查通过；
- 伴随梯度检查通过；
- Gauss--Newton 方向检查通过。

## 10. 兼容性与性能

- 默认选项为 `'auto'`，已有调用代码不需要修改。
- 装配压力均值约束的成本可以忽略。
- 自动判断只扫描一次 `bdFlag`，不增加可见运行成本。
- 不执行矩阵秩检测，因此不会增加大规模问题的求解开销。

## 11. 风险

- `bdFlag==2` 被解释为绝对牵引边界；若调用者使用类型 2 表示其他非标准条件，应显式设置压力约束模式。
- 某些混合边界问题可能存在特殊压力零空间，自动规则不能代替完整的数学适定性分析。
- 显式 `'mean-zero'` 允许用户构造与牵引不相容的问题，这是有意保留的高级覆盖能力。

## 12. 决策记录

1. 保留 `ice_slab` 的顶部零牵引，不通过修改 `g_N` 修复。
2. 修改通用求解器，而不是只修改示例。
3. 默认采用基于牵引边界的自动压力规范。
4. 增加 `'mean-zero'` 和 `'none'` 显式覆盖模式。
5. 自动判断使用所有 `bdFlag==2` 边界。
6. 不使用运行时矩阵秩检测。
7. 正问题和一致切线系统复用同一个约束矩阵。
8. 在 `info` 中公开实际采用的压力约束行为。
