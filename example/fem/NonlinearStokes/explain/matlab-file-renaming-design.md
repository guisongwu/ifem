# Nonlinear Stokes MATLAB 文件重命名设计

## 目标

仅重命名 `example/fem/NonlinearStokes/` 中的 `.m` 文件，采用简洁、概括、
符合 MATLAB/iFEM 习惯的 CamelCase 名称。

## 文件映射

| 旧文件名 | 新文件名 |
|---|---|
| `NonlinearStokesP2P1.m` | `NonlinearStokesP2P1.m` |
| `NSSlab.m` | `NSSlab.m` |
| `NSSlabContinuation.m` | `NSSlabContinuation.m` |
| `NSMMSData.m` | `NSMMSData.m` |
| `NSMMS.m` | `NSConverRate.m` |
| `NSMMSContinuation.m` | `NSEpsContinuation.m` |
| `NSBetaInversion.m` | `NSBetaInversion.m` |
| `NSAdjointInversion.m` | `NSAdjointInversion.m` |
| `NSDiagnosis.m` | `NSDiagnosis.m` |
| `NSRegression.m` | `NSRegression.m` |

## 规则

- 函数文件中的主函数名与文件名同步修改。
- 脚本的首行 section 标题同步修改。
- MATLAB 调用引用和文档正文中的 `.m` 文件引用同步更新。
- 不保留旧文件名包装入口。
- Markdown 与 `.mat` 文件名不修改。
- 数学模型、参数和运行行为不修改。

## 验证

重命名后运行：

```matlab
NSRegression
```

并确认 MATLAB `which` 解析到所有新文件名。

## 决策记录

1. 使用 MATLAB/iFEM CamelCase 风格。
2. 所有文件使用 `NS` 前缀。
3. 周期性不再写入求解器文件名，由选项表达。
4. 不保留旧入口。
5. 本轮不合并或重命名现有文档。
