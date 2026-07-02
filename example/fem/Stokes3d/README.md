# Nonlinear Stokes 3-D

本目录包含三维非线性全 Stokes 冰流模型的 P2--P1 有限元求解器、
制造解收敛测试和 ISMIP-HOM A/C 正问题脚本。

## 快速开始

在 MATLAB 中运行：

```matlab
cd('/path/to/ifem')
setpath
cd example/fem/Stokes3d
```

建议先运行：

```matlab
NSConverRate
ISMIPHOM_A
ISMIPHOM_C
```

## 核心求解器

### `NonlinearStokes3P2P1.m`

通用三维非线性 Stokes 求解器：

```matlab
[soln,eqn,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option);
```

支持：

- P2 速度与 P1 压力 Taylor--Hood 元；
- Glen 非线性黏度；
- Weertman 底部滑移；
- 底面不可穿透约束；
- x、y 方向周期边界；
- 倾斜长方体区域，周期匹配使用 `(x,y,z+sx*x+sy*y)`；
- 可选一致切线矩阵和底部摩擦参数方向导数。

边界标记约定：

- `1`：底面或其他边界上的 no-slip Dirichlet；
- `2`：给定牵引，常用于自由表面；
- `3`：不可穿透 Weertman 滑移底面。

常用数据：

```matlab
pde.A
pde.n
pde.beta
pde.m
pde.f
pde.rho
pde.gravity
pde.g_N
```

常用选项：

```matlab
option.periodic
option.periodic_x
option.periodic_y
option.periodic_slope
option.eps_reg
option.maxIt
option.tol
option.residual_tol
option.damping
option.quadorder
option.facequadorder
option.assemble_tangent
option.pressure_constraint
```

三维周期牵引问题的压力只差一个常数。当前求解器默认
`option.pressure_constraint = 'none'`，不额外添加压力均值乘子约束。

## 制造解验证

### `NSMMSData.m`

给出三维制造解、体力和顶部牵引。该解在 y 方向周期且独立于 y，
可用于检查三维装配、倾斜周期匹配和底面约束。

### `NSConverRate.m`

在多层网格上计算速度和压力的 \(L^2\) 误差并报告收敛阶。压力误差会先
扣除常数均值差。

解析解、区域、边界条件、正则化参数和当前收敛结果见 `NSMMS.md`。

## ISMIP-HOM 示例

### `ISMIPHOM_A.m`

计算 ISMIP-HOM experiment A：倾斜周期区域、正弦床、底面 no slip、
自由表面。

### `ISMIPHOM_C.m`

计算 ISMIP-HOM experiment C：倾斜周期区域、平行上下表面、底部空间变化
摩擦系数、自由表面。

两个脚本默认只跑 \(L=5\) km 的小网格算例；将 `lengthList` 改为
`1000*[5;10;20;40;80;160]` 可做完整长度扫描。脚本按 benchmark 论文的
表 4/5 口径，在自由表面 `y=L/4` 截面统计水平速度范数。
