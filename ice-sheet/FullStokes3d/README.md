# FullStokes3d

本目录包含三维非线性全 Stokes 冰流模型、底部摩擦系数反演、
ISMIP-HOM A/C 基准算例以及制造解验证。主求解器采用 P2--P1
Taylor--Hood 有限元离散，并支持 Glen 非线性黏度、底部滑移、
不可穿透约束和水平周期边界。

## 目录结构

```text
FullStokes3d/
  NonlinearStokes3P2P1.m       % 三维全 Stokes 主求解器
  FullStokes3AdjInvSlab.m      % 无量纲平板床反演
  FullStokes3AdjInvSin.m       % 无量纲正弦床反演
  FullStokes3AdjInvSlabISM.m   % ISM 内部量纲平板床反演
  FullStokes3AdjInvSinISM.m    % ISM 内部量纲正弦床反演
  HOM/                          % ISMIP-HOM A/C 正问题和 C 反演
  test/                         % 制造解、收敛测试及说明
  archive/                      % 已停用或备份脚本，仅供参考
  output/                       % 计算结果和导出的 EPS 图片
```

## 快速开始

```matlab
cd('/path/to/ifem')
setpath
seticepath
cd ice-sheet/FullStokes3d
```

建议先运行：

```matlab
run('test/NSConverRate.m')
run('HOM/HOM_A_L5.m')
run('HOM/HOM_C_L5.m')
```

主要反演入口为：

```matlab
FullStokes3AdjInvSlab
FullStokes3AdjInvSin
FullStokes3AdjInvSlabISM
FullStokes3AdjInvSinISM
```

## 量纲约定

`FullStokes3AdjInvSlab.m` 和 `FullStokes3AdjInvSin.m` 使用无量纲参数，
主要用于验证反演算法。带 `ISM` 后缀的脚本使用 PHGISM/ISM 内部量纲：

- 坐标：km；
- 速度：m/year；
- 底部摩擦系数：Pa yr m^{-1}；
- 压力自由度：物理压力除以 `1e5`。

ISM 脚本通过 `pde.beta_scale` 和 `pde.pressure_dof_scale` 将参数传给
同一个 `NonlinearStokes3P2P1.m` 求解器，因此两类脚本只在参数量纲和
缩放上不同，不使用两套求解器。

## 主求解器

### `NonlinearStokes3P2P1.m`

```matlab
[soln,eqn,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option);
```

主要支持：

- P2 速度与 P1 压力 Taylor--Hood 元；
- Glen 非线性黏度；
- Weertman 底部滑移和底面不可穿透约束；
- x、y 方向周期边界；
- 倾斜长方体区域的周期匹配；
- 一致切线矩阵和底部摩擦参数方向导数。

边界标记约定：

- `1`：no-slip Dirichlet 边界；
- `2`：给定牵引，通常用于自由表面；
- `3`：不可穿透 Weertman 滑移底面。

常用数据字段为 `pde.A`、`pde.n`、`pde.beta`、`pde.m`、`pde.f`、
`pde.rho`、`pde.gravity` 和 `pde.g_N`。常用选项包括周期边界、
正则化黏度、非线性迭代容差、积分阶数、一致切线矩阵和压力约束。

## 反演脚本

`FullStokes3AdjInvSlab.m` 和 `FullStokes3AdjInvSin.m` 分别在平板床与
正弦床区域上进行无量纲伴随反演。对应的 `ISM` 脚本采用冰盖模型内部
量纲，可用于和 PHGISM 的变量及结果直接比较。

所有反演脚本将 EPS 图片写入：

```text
output/<scriptName>/
```

## ISMIP-HOM 基准

`HOM/` 目录包含：

```text
HOM_A.m             % experiment A 多长度正问题
HOM_A_L5.m          % experiment A，L = 5 km 场变量图
HOM_C.m             % experiment C 多长度正问题
HOM_C_L5.m          % experiment C，L = 5 km 场变量图
HOM_CAdjInv.m       % experiment C 无量纲反演
HOM_CAdjInvISM.m    % experiment C ISM 内部量纲反演
```

从 `FullStokes3d/` 中可以直接运行：

```matlab
run('HOM/HOM_A.m')
run('HOM/HOM_C.m')
run('HOM/HOM_CAdjInvISM.m')
```

## 测试与归档

`test/NSMMSData.m` 给出三维制造解、体力和顶部牵引，
`test/NSConverRate.m` 计算速度和压力误差及收敛阶，详细推导见
`test/NSMMS.md`。修改主求解器后至少运行：

```matlab
run('test/NSConverRate.m')
```

`archive/` 中保存早期矩形区域示例、备份脚本和已被当前入口替代的代码。
这些文件不作为测试或推荐入口，也不会自动加入正常计算流程。
