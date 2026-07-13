# 3-D Manufactured-Solution Test

本文档说明 `NSMMSData.m` 和 `NSConverRate.m` 使用的三维解析解、计算区域、
边界条件、正则化参数和收敛阶结果。

## 求解区域

收敛测试使用参考区域

```text
(xi,eta,zeta) in [0,1] x [0,1] x [0,1].
```

物理区域由参考区域映射得到：

```text
x = L xi,
y = W eta,
z = -slope*x + H zeta.
```

当前测试参数为

```text
L = 1,
W = 1,
H = 0.5,
slope = 0.1.
```

因此底面为 `q = 0`，自由表面为 `q = H`，其中

```text
q = z + slope*x.
```

网格由 `cubemesh([0,1,0,1,0,1],h)` 在参考立方体上生成，然后直接映射到
上述斜板区域。当前收敛测试不再扰动内部节点，因此使用规则四面体网格。

## 解析解

令

```text
kx = 2*pi/L,
ky = 2*pi/W,
F(q) = q^2 sin(pi*q/H),
P(q) = 1 + q + 0.3 sin(pi*q/H).
```

定义流函数

```text
psi(x,y,z) = sin(kx*x) sin(ky*y) F(q).
```

速度取为向量势 `(0,psi,0)` 的旋度：

```text
u = curl(0,psi,0) = (-psi_z, 0, psi_x).
```

展开为

```text
u_x = -sin(kx*x) sin(ky*y) F'(q),
u_y = 0,
u_z =  kx cos(kx*x) sin(ky*y) F(q)
     + slope sin(kx*x) sin(ky*y) F'(q).
```

该速度严格满足

```text
div u = 0.
```

压力为

```text
p = cos(kx*x) cos(ky*y) P(q).
```

这里使用 Glen 非线性黏度，参数为：

```text
A = 1,
n = 3.
```

正则化黏度与 `NonlinearStokes3P2P1.m` 中的定义一致：

```text
eta(u) = 0.5 A^(-1/n) (epsilon_II(u) + eps_reg^2)^((1-n)/(2n)).
```

体力 `pde.f` 按

```text
f = -div(2 eta epsilon(u)) + grad p
```

生成。程序中 `eta` 和 `epsilon(u)` 用解析公式评价，`div(2 eta epsilon(u))`
用四阶中心差分对解析应力求导。顶部牵引 `pde.g_N` 按

```text
sigma n,  sigma = 2 eta epsilon(u) - p I
```

精确计算。

## 边界条件

边界标记由 `NSConverRate.m` 设置：

```matlab
bdFlag = setboundary3(node,elem,'Neumann','z==1',...
    'Dirichlet','z==0');
```

对应物理边界条件为：

- `x=0` 与 `x=L`：周期边界；
- `y=0` 与 `y=W`：周期边界；
- 顶面 `q=H`：给定精确牵引 `sigma n = g_N`；
- 底面 `q=0`：no-slip Dirichlet，`u=0`。

由于 `F(0)=0` 且 `F'(0)=0`，解析速度在底面严格为零。

周期匹配使用求解器选项

```matlab
option.periodic = true;
option.periodic_x = [0,L];
option.periodic_y = [0,W];
option.periodic_slope = [slope,0];
```

即在坐标 `(x,y,z+slope*x)` 中匹配周期面。

## 正则化参数

收敛测试传入

```matlab
eps_reg = 1e-2;
option.eps_reg = eps_reg;
option.tol = 1e-8;
option.residual_tol = 1e-8;
```

当前 MMS 使用 `n=3`，因此该正则化参数直接进入 Glen 黏度：

```text
eta(u) = 0.5 (epsilon_II(u) + eps_reg^2)^(-1/3).
```

该参数同时用于制造体力和数值求解器，保证二者一致。

压力常数不固定：

```matlab
option.pressure_constraint = 'none';
```

计算压力误差时，脚本先扣除数值压力与精确压力的体积平均差。

## 收敛阶结果

运行命令：

```matlab
cd('/path/to/ifem')
setpath
seticepath
cd ice-sheet/FullStokes3d
run('test/NSConverRate.m')
```

`NSConverRate.m` 当前默认网格层次为

```matlab
hlist = [1/2;1/4;1/8];
```

由于相邻网格尺寸均为二分加密，脚本用

```matlab
log2(err(k)/err(k+1))
```

计算相邻两层的收敛阶。

其中 `n=3` 后非线性 Picard 迭代明显更慢。当前已经完成的 n=3 检查结果为：

| h       | velocity L2 | u rate | pressure L2 | p rate | Picard it |
| :---:   | :---:       | :---:  | :---:       | :---:  | :---:     |
| 0.50000 | 1.5070e-01  | --     | 3.9420e-01  | --     | 61        |
| 0.25000 | 3.6760e-02  | 2.04   | 1.5970e-01  | 1.30   | 62        |
| 0.12500 | 4.1482e-03  | 3.15   | 3.1674e-02  | 2.33   | 62        |

P2--P1 Taylor--Hood 元的典型期望为

```text
||u-u_h||_L2 = O(h^3),
||p-p_h||_L2 = O(h^2).
```

当前已完成的两层 n=3 检查只能作为规则网格流程的 smoke test。完整定量判断
应以三层或更多加密网格为准。
