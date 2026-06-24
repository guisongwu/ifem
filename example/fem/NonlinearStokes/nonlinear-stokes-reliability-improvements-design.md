# 非线性 Stokes 可靠性改进设计

## 1. 目标

本设计覆盖以下五项改进：

1. 保留并明确记录 `pressure_constraint='auto'` 的边界判断结果；
2. 增加无牵引边界的适定性回归测试；
3. 建立完整 MATLAB 回归测试套件；
4. 验证有限差分反演的稳定性和目标函数下降；
5. 使用完整非线性离散残差参与 Picard 停止判断。

## 2. 范围

主要修改：

```text
example/fem/NonlinearStokes/NonlinearStokesP2P1.m
```

新增完整回归测试：

```text
example/fem/NonlinearStokes/NonlinearStokesRegression.m
```

不修改：

- 旧线性求解器 `StokesP2P1_periodic.m`；
- Glen 黏度模型；
- Weertman 底部滑移模型；
- 周期边界和底部不可穿透条件；
- 反演观测模型。

## 3. 压力自动约束

保留现有规则：

```matlab
hasTractionBoundary = any(bdFlag(:)==2);
```

在 `auto` 模式中：

```text
存在牵引边界 → 不固定压力均值
无牵引边界   → 固定压力零均值
```

不增加运行时秩检测。

继续通过 `info` 报告：

```matlab
info.pressureConstraint
info.hasTractionBoundary
info.pressureMeanConstrained
```

## 4. 完整非线性残差

状态和约束变量记为

$$
X=
\begin{bmatrix}
u\\
p
\end{bmatrix},
\qquad
\lambda=\text{约束拉格朗日乘子}.
$$

在当前状态 $u$ 重新装配：

$$
K(u),\qquad K_b(u),\qquad F.
$$

完整增广离散残差定义为

$$
R_{\mathrm{momentum}}
=
\left(K(u)+K_b(u)\right)u
+B^Tp
-F
+C_u^T\lambda,
$$

$$
R_{\mathrm{div}}
=Bu+C_p^T\lambda,
$$

$$
R_{\mathrm{constraint}}
=CX.
$$

其中将约束矩阵按速度和压力列分块：

$$
C=[C_u,\ C_p].
$$

也可统一写为

$$
R_{\mathrm{state}}
=
\begin{bmatrix}
K(u)+K_b(u) & B^T\\
B & 0
\end{bmatrix}
X
-
\begin{bmatrix}
F\\
0
\end{bmatrix}
+C^T\lambda.
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

必须保留 $\lambda$，否则周期约束、底部不可穿透约束和压力规范产生的反力会被错误计入动量或连续性残差。

## 5. 阻尼与拉格朗日乘子

当前 Picard 迭代对速度和压力使用阻尼：

$$
u^{k+1}=(1-\alpha)u^k+\alpha\widehat u^{k+1},
$$

$$
p^{k+1}=(1-\alpha)p^k+\alpha\widehat p^{k+1}.
$$

约束乘子也应同步更新：

$$
\lambda^{k+1}
=(1-\alpha)\lambda^k+\alpha\widehat\lambda^{k+1}.
$$

这样残差检查对应同一个阻尼后的增广状态。

初始乘子取零。

## 6. 残差归一化

定义状态右端：

$$
b=
\begin{bmatrix}
F\\
0
\end{bmatrix}.
$$

状态残差归一化为

$$
r_{\mathrm{state}}
=
\frac{\|R_{\mathrm{state}}\|}
{\max(1,\|b\|)}.
$$

约束残差归一化为

$$
r_{\mathrm{constraint}}
=
\frac{\|CX\|}
{\max(1,\|X\|)}.
$$

总残差定义为

$$
r_{\mathrm{total}}
=\max(r_{\mathrm{state}},r_{\mathrm{constraint}}).
$$

同时单独记录：

- 动量残差；
- 连续性残差；
- 约束残差；
- 总残差。

## 7. 新增选项

新增：

```matlab
option.residual_tol
```

默认：

```matlab
option.residual_tol = option.tol;
```

新增：

```matlab
option.residual_check_threshold
```

默认：

```matlab
option.residual_check_threshold = ...
    max(1e-2,sqrt(option.residual_tol));
```

该阈值控制何时开始每步装配完整残差。

## 8. 延迟残差检查

为了避免每次 Picard 迭代都额外装配非线性矩阵：

1. 每步计算相对增量；
2. 当

$$
\mathrm{relativeChange}
\le
\mathrm{residualCheckThreshold}
$$

时开始完整残差检查；
3. 一旦开始，之后每步都检查；
4. 若达到最大迭代数，最后一步必须检查；
5. 只有同时满足

$$
\mathrm{relativeChange}<\mathrm{tol}
$$

和

$$
r_{\mathrm{total}}<\mathrm{residual\_tol}
$$

才判定收敛。

未检查的历史位置记录为 `NaN`。

## 9. `info` 输出

新增：

```matlab
info.momentumResidual
info.divergenceResidual
info.constraintResidual
info.nonlinearResidual
info.residualChecked
info.residualTolerance
```

其中每个残差历史长度与 `info.relchange` 一致。

最终判定：

```matlab
info.converged = ...
    info.relchange(end) < option.tol && ...
    info.nonlinearResidual(end) < option.residual_tol;
```

## 10. 无牵引边界测试

构造周期关闭或具有充分速度边界约束的 Stokes 测试，使：

```matlab
any(bdFlag(:)==2) == false
```

验证：

- `auto` 选择压力零均值；
- `info.pressureMeanConstrained == true`；
- 压力均值接近零；
- 增广矩阵满秩；
- 非线性残差达到容差。

测试不追求复杂物理解，只验证压力常数零空间处理。

## 11. 完整回归测试

新增：

```text
example/fem/NonlinearStokes/NonlinearStokesRegression.m
```

允许运行数分钟，覆盖以下内容。

### 11.1 `ice_slab` 三层网格

使用：

```matlab
h = 1/8, 1/16, 1/32
```

检查：

- 正问题收敛；
- `auto` 不固定压力均值；
- 散度积分接近零；
- 散度 $L^2$ 范数随网格下降；
- 完整非线性残差达到容差；
- 周期和底部约束残差达到容差。

### 11.2 MMS 收敛阶

检查：

- 速度 $L^2$ 收敛阶不低于合理阈值；
- 压力 $L^2$ 收敛阶不低于合理阈值；
- 散度误差随网格下降；
- 每层完整残差达到容差。

### 11.3 压力模式

检查：

- `ice_slab` 中 `auto` 与 `none` 一致；
- `mean-zero` 明确启用均值约束；
- MMS 中 `auto` 与 `mean-zero` 一致；
- 未知选项产生预期错误。

### 11.4 直接与延拓

检查 `ice_slab`：

- 直接求解收敛；
- 延拓各阶段收敛；
- 最终速度和压力在容差内一致。

### 11.5 无牵引边界

执行第 10 节的适定性测试。

### 11.6 伴随反演

检查：

- 状态导数相对误差；
- 伴随梯度相对误差；
- Gauss--Newton 方向误差；
- PCG 正常收敛；
- 目标函数下降。

为便于自动测试，反演脚本中的导数检查应返回数值，而不仅输出文本。

### 11.7 有限差分反演

稳定性验收标准：

- 所有正问题收敛；
- 目标函数为有限实数；
- 接受更新后的目标函数不增加；
- 参数保持正值；
- 迭代历史中无 `NaN` 或 `Inf`。

不要求最终参数误差达到指定高精度。

## 12. 性能

- 残差检查开始前不增加额外非线性装配。
- 接近收敛后每步增加一次黏性、底部和载荷装配。
- 完整回归测试允许运行数分钟。
- 不进行运行时矩阵秩检测。

## 13. 可靠性

- 正问题只有在增量和完整残差都满足容差时才报告收敛。
- 回归测试任一断言失败即停止并给出明确错误。
- 反演测试以稳定性和下降性为核心，不把病态反问题的恢复精度误判为求解器正确性。

## 14. 维护

- 测试脚本使用局部辅助函数，避免污染库 API。
- 关键容差集中定义，便于未来调整。
- 回归脚本应输出紧凑摘要，详细迭代输出默认关闭。
- MATLAB R2025b 是本轮主要验证环境。
- 本轮不要求 Octave 兼容。

## 15. 决策记录

1. `pressure_constraint='auto'` 保持现有边界标记规则。
2. 不增加数值秩检测或自动适定性猜测。
3. 新增无牵引边界回归测试。
4. 建立允许数分钟运行的完整 MATLAB 回归套件。
5. 有限差分反演只要求稳定、有限和目标下降。
6. 非线性收敛必须同时满足增量与完整残差。
7. 完整残差包含拉格朗日乘子反力。
8. 拉格朗日乘子与速度、压力同步阻尼更新。
9. 使用延迟残差装配降低运行成本。
10. 不修改旧线性 Stokes 求解器。
