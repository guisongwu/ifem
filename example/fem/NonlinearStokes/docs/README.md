# Nonlinear Stokes 说明文档

本目录按用途整理非线性 Stokes 示例的说明文档：

- `theory/`：连续数学模型、伴随方程、压力规范等理论说明；
- `implementation/`：求解器、切线矩阵、残差检查和反演脚本的实现说明；
- `design/`：迁移、重命名、延拓、压力约束和回归测试等设计记录。

## Theory

- [equations.md](theory/equations.md)：线性滑移情形下正问题、增量正问题、伴随问题和增量伴随问题的连续形式。
- [stress-linearization.md](theory/stress-linearization.md)：通用滑移指数 `m` 下的应力线性化和反问题方程。
- [pressure-traction.md](theory/pressure-traction.md)：压力零均值规范与牵引边界的关系。

## Implementation

- [tangent-matrix.md](implementation/tangent-matrix.md)：残差、黏性残差、一致切线矩阵和 Picard 矩阵的区别。
- [residual-damping-gradient.md](implementation/residual-damping-gradient.md)：严格残差检查、乘子阻尼和伴随梯度有限差分检查。
- [boundary-objective.md](implementation/boundary-objective.md)：边界观测目标函数反演脚本详解。

## Design

- [fd-newton.md](design/fd-newton.md)：有限差分 Newton 反演设计。
- [continuation.md](design/continuation.md)：正则化参数延拓示例设计。
- [migration.md](design/migration.md)：从 `Stokes/` 迁移到独立目录的设计。
- [pressure-constraint.md](design/pressure-constraint.md)：压力约束自动判断与修复设计。
- [reliability.md](design/reliability.md)：可靠性和回归测试改进设计。
- [renaming.md](design/renaming.md)：MATLAB 文件重命名设计。
