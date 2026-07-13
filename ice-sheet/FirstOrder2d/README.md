# FirstOrder2d

本目录包含二维 x-z 截面的 first-order/Blatter--Pattyn 冰流求解器和
底部摩擦系数反演脚本。速度采用 P2 三角形有限元离散，支持 Glen
非线性黏度、周期侧边界、自由表面牵引和底部滑移。

## 目录结构

```text
FirstOrder2d/
  FirstOrderP2.m              % 二维 FO 主求解器
  FirstOrderAdjInvSlab.m      % 平板床伴随反演
  FirstOrderAdjInvSlabISM.m   % PHGISM 单位平板床伴随反演
  FirstOrderAdjInvSin.m       % 正弦床伴随反演
  FirstOrderAdjInvSinISM.m    % PHGISM 单位正弦床伴随反演
  test/                        % 收敛率、可视化和反演对照实验
  docs/                        % 方程和离散说明
  output/                      % 生成的图片和已有计算结果
  README.md
```

## 快速开始

```matlab
cd('/path/to/ifem')
setpath
seticepath
cd ice-sheet/FirstOrder2d
```

推荐入口：

```matlab
FirstOrderAdjInvSlab
FirstOrderAdjInvSlabISM
FirstOrderAdjInvSin
FirstOrderAdjInvSinISM
run('test/FirstOrderConverRate.m')
```

主求解器接口为：

```matlab
[soln,eqn,info,node,elem,bdFlag] = FirstOrderP2(option);
```

## 反演和测试

`FirstOrderAdjInvSlab.m` 使用周期倾斜平板区域，
`FirstOrderAdjInvSin.m` 使用类似 ISMIP-HOM-B 的正弦床区域。
`FirstOrderAdjInvSlabISM.m` 和 `FirstOrderAdjInvSinISM.m` 分别是对应的
PHGISM 单位量纲版本：坐标为 km，速度为 m/year，反演变量为物理底部
摩擦系数 `\beta`。这些脚本均由顶部水平速度观测反演底部摩擦系数。

`test/` 中包含：

```text
FirstOrderConverRate.m              % 制造解收敛率
FirstOrderVisualize.m               % 二维场可视化
FirstOrderAdjInvSinNoRegular.m      % 无正则化反演对照
FirstOrderAdjInvSinParamSweep.m     % 正则化和算法参数扫描
```

所有反演脚本将 EPS 图片写入 `output/<scriptName>/`。迁移前已经生成的
EPS 结果也保留在 `output/` 中。方程、伴随方法和离散细节见
`docs/equations.md`。

## 结构决策

主目录仅保留稳定求解器和正式反演入口；验证、可视化和参数研究归入
`test/`；文档与生成结果分别归入 `docs/` 和 `output/`。此次重构不改变
模型量纲、数值参数、网格或求解算法。
