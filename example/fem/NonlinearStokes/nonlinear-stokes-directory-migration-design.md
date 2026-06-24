# Nonlinear Stokes 独立目录迁移设计

## 目标

将非线性 Stokes 求解器、示例、反演、诊断、测试和文档从

原目录为 `example/fem/Stokes/`。

迁移到与 `Stokes/` 并列的独立目录：

```text
example/fem/NonlinearStokes/
```

## 迁移范围

迁移：

- `NonlinearStokesP2P1.m`
- `NonlinearStokesSlab.m`
- `NonlinearStokesSlabContinuation.m`
- `NonlinearStokesMMSData.m`
- `NonlinearStokesMMS.m`
- `NonlinearStokesMMSContinuation.m`
- `NonlinearStokesBetaInversion.m`
- `NonlinearStokesAdjointInversion.m`
- `NonlinearStokesDiagnosis.m`
- `diagnose_nonlinear_stokes_inversion_result.mat`
- `NonlinearStokesRegression.m`
- 非线性 Stokes 相关 Markdown 文档。

保留在 `Stokes/`：

- 旧线性 Stokes 求解器；
- 旧线性反演脚本；
- `codesinmypaper-explained.md`；
- 其他线性 Stokes 数据、图像和调试文件。

## 目录结构

第一阶段保持扁平结构，不再划分 `solver/`、`inverse/`、`tests/` 或
`docs/` 子目录，以避免增加 MATLAB 路径和文档维护复杂度。

## 兼容性

- MATLAB 文件名和函数名保持不变。
- 脚本入口名称保持不变。
- `setpath` 使用 `genpath`，新目录会自动加入 MATLAB 搜索路径。
- 迁移后更新仓库内所有旧路径引用。

## 验证

迁移完成后执行：

```matlab
NonlinearStokesRegression
```

并检查：

- 正问题；
- MMS；
- 压力约束模式；
- 延拓；
- 伴随导数；
- 有限差分反演。

## 风险

- 若新旧目录存在同名函数，MATLAB 路径顺序可能导致遮蔽。
- Markdown 中的旧路径可能失效。
- 脚本若依赖当前工作目录，迁移后可能暴露隐含相对路径依赖。

## 决策记录

1. 使用 `example/fem/NonlinearStokes/`，与 `Stokes/` 并列。
2. 只迁移非线性内容。
3. 保持扁平目录。
4. 不重命名 MATLAB 入口。
5. 不迁移解释旧线性代码的文档。
6. 不保留编辑器 `.swp` 临时文件。
7. 完整 MATLAB 回归是迁移验收条件。
