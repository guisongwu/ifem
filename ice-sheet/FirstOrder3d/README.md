# FirstOrder3d

本目录包含三维 first-order/Blatter--Pattyn 冰流模型、P2 四面体有限元
求解器和底部摩擦系数伴随反演。模型求解两个水平速度分量，并支持 Glen
非线性黏度、底部滑移和水平周期边界。

## 目录结构

```text
FirstOrder3d/
  NonlinearFOP2.m             % 三维 FO 主求解器
  FirstOrder3AdjInvSlab.m     % 平板床伴随反演
  FirstOrder3AdjInvSlabISM.m  % PHGISM 单位平板床伴随反演
  FirstOrder3AdjInvSin.m      % 正弦床伴随反演
  FirstOrder3AdjInvSinISM.m   % PHGISM 单位正弦床伴随反演
  test/                        % 收敛率、可视化和周期区域诊断
  output/                      % 生成的图片和计算结果
  README.md
```

## 快速开始

```matlab
cd('/path/to/ifem')
setpath
seticepath
cd ice-sheet/FirstOrder3d
```

推荐入口：

```matlab
FirstOrder3AdjInvSlab
FirstOrder3AdjInvSlabISM
FirstOrder3AdjInvSin
FirstOrder3AdjInvSinISM
run('test/FirstOrder3ConverRate.m')
run('test/FirstOrder3ConverRateN3.m')
```

主求解器接口为：

```matlab
[soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option);
```

`pde` 提供 Glen 参数、体力、边界速度、顶部牵引和底部摩擦系数；
`option` 控制周期边界、正则化黏度、非线性迭代和积分阶数。

## 反演和测试

`FirstOrder3AdjInvSlab.m` 和 `FirstOrder3AdjInvSin.m` 分别在三维平板床
与正弦床区域上，从顶部水平速度观测反演底部摩擦系数。
`FirstOrder3AdjInvSlabISM.m` 和 `FirstOrder3AdjInvSinISM.m` 是对应的
PHGISM 单位量纲版本：坐标为 km，速度为 m/year，反演变量为物理底部
摩擦系数 `\beta`。

`test/` 中包含：

```text
FirstOrder3ConverRate.m          % n = 1 制造解收敛率
FirstOrder3ConverRateN3.m        % n = 3 非线性制造解收敛率
FirstOrder3Visualize.m           % 制造解场可视化
FirstOrder3Periodic.m            % 周期长方体诊断求解器
FirstOrder3PeriodicVisualize.m   % 周期长方体结果可视化
```

反演脚本将 EPS 图片写入 `output/<scriptName>/`。可视化脚本启用导出选项后，
也默认写入 `output/<scriptName>/`。

## 结构决策

主目录仅保留通用求解器和正式反演入口；制造解、可视化及特定区域诊断归入
`test/`；所有生成结果集中到 `output/`。当前没有独立的 ISMIP-HOM 数据或
归档代码，因此不创建空的 `HOM/` 或 `archive/`。此次重构不改变模型量纲、
数值参数、网格或求解算法。
