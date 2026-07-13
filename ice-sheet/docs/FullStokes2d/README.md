# Nonlinear Stokes 说明文档

本目录整理 `ice-sheet/FullStokes2d/` 中非线性 Stokes 正问题、伴随方程和底部滑移反演脚本的说明文档。当前文档分为两类：

- 根目录文档：连续模型、符号约定、伴随推导和压力约束；
- `implementation/`：MATLAB 实现、离散切线矩阵、残差检查和反演脚本说明。

## 建议阅读顺序

1. 先读 [tensor-notation.md](tensor-notation.md)，了解张量、切向投影和四阶切线张量记号。
2. 再读 [equations.md](equations.md)，建立正问题、增量正问题、伴随问题和增量伴随问题的连续形式。
3. 如果关注通用 Weertman 滑移指数，读 [weertman-linearization.md](weertman-linearization.md)。
4. 如果关注代码实现和反演脚本，读 `implementation/` 下的文档。

## 连续模型与推导

- [tensor-notation.md](tensor-notation.md)：非线性 Stokes 推导中使用的张量积、双点积、四阶张量和转置记号。
- [equations.md](equations.md)：线性滑移 `m=1` 情形下正问题、增量正问题、伴随问题和增量伴随问题的连续形式。
- [adjoint.md](adjoint.md)：从增量正问题推导伴随问题、伴随边界条件和梯度公式。
- [incremental-adjoint.md](incremental-adjoint.md)：Gauss--Newton Hessian-vector product 中增量伴随问题的来源。
- [weertman-linearization.md](weertman-linearization.md)：通用滑移指数 `m` 下底部 Weertman 滑移项的线性化。
- [pressure-constrain.md](pressure-constrain.md)：压力零均值规范、牵引边界和求解器压力约束之间的关系。

## MATLAB 实现说明

- [implementation/tangent-matrix.md](implementation/tangent-matrix.md)：离散残差、黏性残差、一致切线矩阵和 Picard 矩阵的区别。
- [implementation/residual-damping-gradient.md](implementation/residual-damping-gradient.md)：严格残差检查、乘子阻尼和伴随梯度有限差分检查。
- [implementation/boundary-objective.md](implementation/boundary-objective.md)：`NonlinearStokesAdjInvSlabBed.m` 的边界观测目标函数、上表面线积分、伴随梯度和矩阵自由 Gauss--Newton 实现。
