---
title: codesinmypaper.m：线性 Stokes 底部参数反演逐行解析
tags:
  - 冰盖反演
  - Stokes 方程
  - 伴随方法
  - MATLAB
createTime: 2026/06/22 16:00:00
permalink: /posts/codesinmypaper-linear-stokes-inversion/
---

# `codesinmypaper.m`：线性 Stokes 底部参数反演逐行解析

本文逐行解释 iFEM 示例目录中的 `codesinmypaper.m`。这个脚本在一个倾斜、左右周期的二维冰层截面上求解线性 Stokes 方程，根据上表面的水平速度观测反演底部 Robin 滑移系数。

脚本的主线是

$$
m
\longrightarrow
(\boldsymbol u,p)
\longrightarrow
u_x|_{\Gamma_t}
\longrightarrow
\Xi(m),
$$

然后通过伴随方程计算梯度，并通过 Gauss--Newton Hessian 对参数增量的作用，使用 Krylov 迭代求出参数更新量。

需要先说明：代码中的参数名 `m` 实际表示 Robin 系数或底部摩擦系数，数学上更常记作 $\beta$。它不是 Glen 流动定律中的指数，也不是 Weertman 定律中的幂指数。

## 1. 计算区域

程序先在参考矩形

$$
\widehat\Omega=[0,1]\times[0,0.5]
$$

上生成三角网格，然后执行坐标变换

$$
y\leftarrow y-sx,
\qquad s=0.1.
$$

因此实际计算区域是平行四边形

$$
\Omega
=
\left\{
(x,y):
0<x<1,\quad -sx<y<0.5-sx
\right\}.
$$

四条边分别为

$$
\Gamma_b=\{(x,y):y=-sx\},
$$

$$
\Gamma_t=\{(x,y):y=0.5-sx\},
$$

以及左右边界 $\Gamma_l$ 和 $\Gamma_r$。左右边界满足周期条件。由于区域倾斜，周期对应点之间存在竖直位移，但在参考坐标

$$
\widehat y=y+sx
$$

中，它们具有相同的 $\widehat y$。

底边和顶边的单位法向量可以分别写成

$$
\boldsymbol n_b
=
\frac{(-s,-1)^T}{\sqrt{1+s^2}},
\qquad
\boldsymbol n_t
=
\frac{(s,1)^T}{\sqrt{1+s^2}}.
$$

底边单位切向量可取

$$
\boldsymbol t
=
\frac{(1,-s)^T}{\sqrt{1+s^2}}.
$$

## 2. 代码求解的线性 Stokes 方程

未知量为二维速度

$$
\boldsymbol u=(u_x,u_y)^T
$$

和压力 $p$。`StokesP2P1_periodic.m` 中装配的是分量 Laplace 型线性 Stokes 系统：

$$
-\Delta\boldsymbol u+\nabla p=\boldsymbol f
\qquad\text{in }\Omega,
$$

$$
\nabla\cdot\boldsymbol u=0
\qquad\text{in }\Omega.
$$

在本算例中，体力来自 `stokes_data_grav_period.m`：

$$
\boldsymbol f=(0,-1)^T.
$$

这可以理解为无量纲化后的重力。

### 2.1 顶部边界

顶部被标记为 Neumann 边界。抽象写法为

$$
\partial_n\boldsymbol u-p\boldsymbol n
=
\boldsymbol g_N
\qquad\text{on }\Gamma_t.
$$

当前数据文件令

$$
\boldsymbol g_N=\boldsymbol 0,
$$

即顶部采用齐次自然边界。

### 2.2 底部边界

底部被标记为 Robin 边界，并且 `option.use_slip=true`。求解器会旋转速度自由度，将其分成法向和切向分量。

底部法向速度满足

$$
\boldsymbol u\cdot\boldsymbol n_b=g_{Dn}.
$$

当前算例中

$$
g_{Dn}=0,
$$

即底部不可穿透。

切向分量满足 Robin 型滑移条件。概念上可写成

$$
\partial_n u_t+m(x)u_t=g_{RN,t}
\qquad\text{on }\Gamma_b,
$$

其中

$$
u_t=\boldsymbol u\cdot\boldsymbol t.
$$

代码把 `pde.g_R` 作为 Robin 系数 $m$，把 `pde.g_RN` 作为 Robin 边界右端。当前数据文件给出零右端，因此核心关系是

$$
\partial_n u_t+m(x)u_t=0.
$$

严格地说，求解器以旋转后的离散速度系统实现滑移约束；上式是理解该实现最直接的连续形式。

### 2.3 周期边界

左右边界满足

$$
\boldsymbol u|_{\Gamma_l}
=
\boldsymbol u|_{\Gamma_r},
\qquad
p|_{\Gamma_l}
=
p|_{\Gamma_r},
$$

这里的“相等”是指参考坐标高度 $\widehat y=y+sx$ 相同的对应点。求解器通过合并周期自由度实现这一条件。

## 3. 弱形式与有限元离散

忽略边界数据的细节，线性 Stokes 弱形式为：寻找
$(\boldsymbol u,p)$，使得对所有测试函数
$(\boldsymbol v,q)$ 有

$$
\int_\Omega
\nabla\boldsymbol u:\nabla\boldsymbol v\,dx
-
\int_\Omega
p\,\nabla\cdot\boldsymbol v\,dx
+
\int_{\Gamma_b}
m\,u_t v_t\,ds
=
\int_\Omega
\boldsymbol f\cdot\boldsymbol v\,dx
+
\int_{\Gamma_t}
\boldsymbol g_N\cdot\boldsymbol v\,ds
+
\int_{\Gamma_b}
g_{RN,t}v_t\,ds,
$$

$$
-\int_\Omega q\,\nabla\cdot\boldsymbol u\,dx=0.
$$

离散空间采用 Taylor--Hood 元：

- 速度使用连续分片二次元 $P_2$；
- 压力使用连续分片一次元 $P_1$；
- 底部参数的独立自由度使用周期 $P_1$ 表示，再插值到 $P_2$ 边中点。

离散线性系统具有鞍点结构

$$
\begin{pmatrix}
A(m) & B^T\\
B & 0
\end{pmatrix}
\begin{pmatrix}
\boldsymbol U\\
\boldsymbol P
\end{pmatrix}
=
\begin{pmatrix}
\boldsymbol F\\
0
\end{pmatrix},
$$

其中

$$
A(m)=A_\Omega+A_{\Gamma_b}(m).
$$

$A_\Omega$ 是速度 Laplace 矩阵，$A_{\Gamma_b}(m)$ 是底部 Robin 矩阵。

## 4. 反问题、目标泛函与观测

脚本先指定“真实参数”

$$
m_0(x)
=
1+0.1\cos(2\pi x+0.1\pi),
$$

求解正问题并把顶部速度保存为合成观测。

虽然 `u_obs` 同时保存了水平和竖直速度，但目标泛函实际只使用顶部水平速度：

$$
\Xi(m)
=
\int_{\Gamma_t}
\left(
u_x(m)-u_x^{\mathrm{obs}}
\right)^2\,ds.
$$

因此

$$
\frac{\partial\Xi}{\partial u_x}
=
2\left(u_x-u_x^{\mathrm{obs}}\right),
\qquad
\frac{\partial\Xi}{\partial u_y}=0.
$$

`integral_neumann_P2(...,'repeat',...)` 用 P2 边界求积计算这个平方误差。

## 5. 伴随梯度

把离散或连续状态方程抽象写成

$$
L(\boldsymbol u,p;m)=0.
$$

参数方向 $\delta m$ 引起的速度增量 $\delta\boldsymbol u$ 满足

$$
L_u\,\delta\boldsymbol u
=
-L_m\,\delta m.
$$

直接逐个参数方向计算 $\delta\boldsymbol u$ 的成本较高。于是引入伴随速度 $\boldsymbol u^*$，使其满足

$$
L_u^T\boldsymbol u^*
=
-\Xi_u^T.
$$

线性 Stokes 算子在当前离散和边界设置下基本对称，因此伴随问题仍调用同一个 Stokes 求解器，只是把顶部边界右端换成

$$
-2
\begin{pmatrix}
u_x-u_x^{\mathrm{obs}}\\
0
\end{pmatrix}.
$$

底部 Robin 项对参数的方向导数为

$$
\int_{\Gamma_b}
\delta m\,
\boldsymbol u\cdot\boldsymbol v\,ds,
$$

或在旋转后的滑移表述中理解为切向速度乘积。于是梯度方向作用为

$$
\Xi'(m)[\delta m]
=
\int_{\Gamma_b}
\delta m\,
\boldsymbol u^*\cdot\boldsymbol u\,ds.
$$

若 $\{\psi_i\}_{i=1}^{N_m}$ 是底部参数基函数，则梯度分量为

$$
g_i
=
\int_{\Gamma_b}
\psi_i\,
\boldsymbol u^*\cdot\boldsymbol u\,ds.
$$

这正是 `integral_robin_P2` 在第 192--198 行中计算的量。

## 6. Gauss--Newton Hessian-vector product

脚本没有显式组装完整 Hessian，而是调用 `stokes_hessian.m` 计算

$$
\boldsymbol h
\longmapsto
H_{\mathrm{GN}}(m)\boldsymbol h.
$$

给定参数方向 $\boldsymbol h$，辅助函数执行两次线性 Stokes 求解。

第一步，求状态增量 $\delta\boldsymbol u$：

$$
L_u\,\delta\boldsymbol u
=
-L_m\boldsymbol h.
$$

第二步，求增量伴随变量 $\boldsymbol w$：

$$
L_u^T\boldsymbol w
=
-2C^TC\,\delta\boldsymbol u,
$$

其中 $C$ 是提取顶部水平速度的观测算子。

最后计算

$$
\left(H_{\mathrm{GN}}\boldsymbol h\right)_i
=
\int_{\Gamma_b}
\psi_i\,
\boldsymbol w\cdot\boldsymbol u\,ds
+
\gamma
\left(M_{\mathrm{stab}}\boldsymbol h\right)_i.
$$

主程序使用 `cgs` 解

$$
H_{\mathrm{GN}}(m)\,\delta m
=
\nabla\Xi(m),
$$

然后更新

$$
m_{\mathrm{new}}
=
m-\delta m.
$$

## 7. MATLAB 数据布局

理解代码前必须先理解三个自由度布局。

### 7.1 P2 速度自由度

若顶点数为 $N$、边数为 $N_E$，则每个速度分量有

$$
N_u=N+N_E
$$

个 P2 自由度。前 $N$ 个对应网格顶点，后 $N_E$ 个对应边中点。

完整速度向量按分量存储：

$$
\texttt{soln.u}
=
\begin{bmatrix}
u_x\\
u_y
\end{bmatrix},
$$

因此

- `1:Nu` 是水平速度；
- `Nu+1:2*Nu` 是竖直速度。

### 7.2 压力自由度

压力使用 P1 元，只定义在顶点上，因此

$$
N_p=N.
$$

### 7.3 底部参数自由度

独立参数只有周期网格顶点上的 $N_m=n_x-1$ 个值。为了供 P2 边界求积和 Robin 矩阵使用，代码把它扩展为

$$
\begin{bmatrix}
\text{底部顶点值}\\
\text{底部边中点值}
\end{bmatrix}.
$$

边中点值由相邻顶点平均得到。

## 8. 第 1--32 行：标题、全局变量和控制参数

### 第 1 行

```matlab
%% Stokes inversion for m
```

`%%` 是 MATLAB 的代码节标记。MATLAB 编辑器可以按节运行脚本。这里说明本节用于反演参数 `m`。

### 第 2 行

空行，只用于提高可读性。

### 第 3 行

```matlab
fprintf('==============Stokes inversion robin beta===============\n');
```

`fprintf` 向命令窗口输出格式化文本。`\n` 表示换行。字符串把反演量称为 Robin 系数 `beta`，进一步说明变量 `m` 就是 $\beta$。

### 第 4 行

```matlab
global slope h dbg_case dbg_on;
```

声明四个全局变量。其他函数只有再次写出同名 `global` 声明后，才能访问这些变量。

- `slope`：倾斜率 $s$；
- `h`：网格尺寸；
- `dbg_case`、`dbg_on`：调试开关。

一般不建议大量使用全局变量，因为它会隐藏函数的数据依赖；这里是旧示例代码的实现方式。

### 第 5 行

空行。

### 第 6--7 行

```matlab
dbg_on = false;
dbg_case = 0;
```

关闭调试输出，并把调试场景编号设为零。`false` 是 MATLAB 逻辑假值。

### 第 8 行

空行。

### 第 9 行

```matlab
%% Control parameters
```

开始“控制参数”代码节。

### 第 10--13 行

```matlab
slope = 0.1;
h = 0.1;
scheme = 6;
max_iteration = 5;
```

- `slope=0.1` 设置区域坡度；
- `h=0.1` 设置均匀网格尺寸；
- `scheme=6` 原本用于在完整版 `stokes_inversion.m` 中选择算法；
- `max_iteration=5` 指定最多执行五次反演更新。

在当前精简脚本中，`scheme` 后面没有被读取，因此它是遗留变量，不会改变算法。

### 第 14 行

空行。

### 第 15--17 行

```matlab
% plot dim
plot_m = 2;
plot_n = 3;
```

`%` 开始单行注释。这里把第一个图窗划分为 $2\times3$ 个子图。

### 第 18--23 行

```matlab
option.verb = 0;
option.solver = 'direct';
option.quadorder = 4;
option.use_newton = false;
option.use_slip = true;
option.periodic = true;
```

MATLAB 结构体可以通过点号动态添加字段。

- `verb=0`：不输出详细求解信息；
- `solver='direct'`：使用反斜杠直接求解离散线性系统；
- `quadorder=4`：区域积分采用指定阶数的三角形求积；
- `use_newton=false`：在该线性求解器路径中没有实际作用；
- `use_slip=true`：启用底部法向/切向旋转和滑移边界处理；
- `periodic=true`：启用左右周期自由度约束。

### 第 24--25 行

空行。

### 第 26 行

```matlab
%% Constants
```

开始常数定义代码节。

### 第 27--30 行

```matlab
four = mp('4');
three = mp('3');
two = mp('2');
one = mp('1');
```

`mp` 通常来自 Multiprecision Computing Toolbox，并不是 MATLAB 基础函数。把字符串转换为多精度数可以减少某些差分和 Hessian 验证中的舍入误差。

当前脚本实际只使用了 `two`；`four`、`three`、`one` 是从实验代码保留下来的变量。

### 第 31--32 行

空行。

## 9. 第 33--61 行：生成并倾斜计算区域

### 第 33--34 行

```matlab
%% Setup domain
figure(1);
```

开始区域设置，并激活编号为 1 的图窗。

### 第 35 行

```matlab
[node,elem] = squaremesh([0 1 0 0.5], h);
```

`squaremesh` 在矩形 $[0,1]\times[0,0.5]$ 上生成尺寸约为 `h` 的三角网格。

MATLAB 的多输出语法

```matlab
[a,b] = function(...)
```

用于同时接收多个返回值：

- `node`：节点坐标矩阵，每行是一个二维坐标；
- `elem`：三角形连接矩阵，每行存放一个三角形的三个顶点编号。

### 第 36 行

```matlab
bdFlag = setboundary(node,elem, ...
    'Neumann','y==0.5', 'Robin', 'y==0');
```

在坐标变换之前标记边界：

- 参考区域的 `y==0.5` 标记为 Neumann；
- 参考区域的 `y==0` 标记为 Robin。

必须先标记再倾斜节点，否则字符串条件 `y==0.5` 和 `y==0` 不再成立。

### 第 37--39 行

```matlab
nx = length([0:h:1]);
ny = length([0:h:.5]);
n1 = nx-1;
```

冒号表达式 `a:h:b` 生成从 `a` 到 `b`、步长为 `h` 的行向量。

- `nx` 是水平方向节点数；
- `ny` 是竖直方向节点数；
- `n1=nx-1` 是周期底边上的独立顶点数。

因为 $x=0$ 和 $x=1$ 周期等价，所以不能把两个端点都当成独立参数自由度。

### 第 40 行

空行。

### 第 41--42 行

```matlab
X = reshape(node(:,1), ny, nx);
Y = reshape(node(:,2), ny, nx);
```

- `node(:,1)` 取 `node` 的全部行、第 1 列；
- `node(:,2)` 取全部竖直坐标；
- `reshape` 在不改变元素总数的条件下重排数组形状。

这里依赖 `squaremesh` 的节点编号顺序，把一维节点列表恢复为 `ny × nx` 的逻辑网格。

### 第 43--44 行

```matlab
assert(norm(X(1,:)' - [0:h:1]') < 1e-12);
assert(norm(Y(:,1) - [0:h:.5]') < 1e-12);
```

`assert(condition)` 在条件为假时终止程序。这里检查节点编号顺序是否符合后续索引假设。

- `X(1,:)` 是第一行；
- `'` 是共轭转置，对实数向量就是转置；
- `norm` 默认计算向量的 Euclidean 范数。

### 第 45 行

空行。

### 第 46--52 行

```matlab
N = size(node, 1);
[elem2dof,edge,bdDof] = dofP2(elem);
NE = size(edge,1);
Nu = N + NE;
Np = N;
Nm = n1;
EI = eye(Nm);
```

- `size(node,1)` 返回矩阵行数，即顶点数；
- `dofP2` 建立 P2 自由度编号；
- `edge` 包含所有网格边；
- P2 标量空间的自由度数是顶点数加边数，所以 `Nu=N+NE`；
- P1 压力只使用顶点，所以 `Np=N`；
- 参数有 `Nm=n1` 个独立周期顶点值；
- `eye(Nm)` 生成 $N_m\times N_m$ 单位矩阵。

`EI(:,ii)` 后面用于取得第 $i$ 个标准基向量。

变量 `elem2dof`、`bdDof` 和 `Np` 在主脚本后续没有直接使用，但前两者有助于说明自由度结构。

### 第 53 行

空行。

### 第 54--61 行

```matlab
if slope ~= 0
    fprintf(2, 'Slab test %f\n', slope);
    node(:,2) = node(:,2) - slope * node(:,1);
end
```

`~=` 表示“不等于”。若坡度非零：

1. `fprintf(2,...)` 向标准错误流输出信息；
2. 把每个节点的第二坐标改为 $y-sx$；
3. `end` 结束 `if` 块。

这一步把矩形变成倾斜冰层区域。

## 10. 第 63--79 行：建立顶部和底部速度索引

### 第 63 行

```matlab
%% Index
```

开始索引设置代码节。

### 第 64--65 行

```matlab
IUxNode = [1:N];
IUyNode = [Nu+1:Nu+N];
```

方括号用于构造数组。由于速度按分量排列：

- `IUxNode` 取得顶点上的 $u_x$；
- `IUyNode` 取得顶点上的 $u_y$。

这里的方括号其实可以省略，`1:N` 本身就是向量。

### 第 66 行

空行。

### 第 67 行

```matlab
unode = [node; ...
    (node(edge(:,1),:)+node(edge(:,2),:))/2];
```

分号 `;` 在方括号内表示竖直拼接。

`unode` 依次存放：

1. 所有 P2 顶点自由度坐标；
2. 所有 P2 边中点自由度坐标。

这与 `dofP2` 的自由度编号顺序一致。

### 第 68 行

```matlab
unode(:,2) = unode(:,2) + slope * unode(:,1);
```

把物理坐标转换回参考高度

$$
\widehat y=y+sx.
$$

这样底边重新对应 $\widehat y=0$，顶边对应 $\widehat y=0.5$，便于筛选自由度。

### 第 69 行

空行。

### 第 70--71 行

```matlab
IUxBot = sort(find(abs(unode(:,2)) < 1e-8 ...
                  & unode(:,1) < 1-h/4));
IUyBot = IUxBot + Nu;
```

- `abs(unode(:,2))<1e-8` 选取参考高度为零的底部自由度；
- `&` 是逐元素逻辑“与”；
- `unode(:,1)<1-h/4` 排除周期右端点 $x=1$；
- `find` 返回逻辑条件为真的索引；
- `sort` 按索引升序排列；
- 竖直速度编号等于水平速度编号加 `Nu`。

使用 `1-h/4` 而不是严格的 `x<1`，是为了在浮点误差下可靠排除右端周期重复点。

### 第 72 行

空行。

### 第 73--74 行

```matlab
IUxTop = sort(find(abs(unode(:,2) - .5) < 1e-8 ...
                  & unode(:,1) < 1-h/4));
IUyTop = IUxTop + Nu;
```

与底部相同，但筛选参考高度 $\widehat y=0.5$ 的顶部 P2 自由度。

### 第 75 行

空行。

### 第 76 行

```matlab
Isft = reshape([1:n1; n1+1:2*n1], 2*n1, 1);
```

若一个周期 P2 边界向量按

```text
[全部顶点值; 全部边中点值]
```

排列，这一索引把它重排为

```text
[顶点1, 中点1, 顶点2, 中点2, ...]
```

以便按几何顺序绘图。

### 第 77--78 行

```matlab
sft = @(dat) dat(Isft);
sft_ext = @(dat) sft(extend_mid(dat));
```

`@(dat)` 创建匿名函数。

- `sft(dat)` 按 `Isft` 重排已有 P2 边界数据；
- `sft_ext(dat)` 先把 P1 参数扩展到 P2，再重排。

`sft_ext` 在当前脚本中没有被调用。

### 第 79 行

空行。

## 11. 第 81--123 行：设置 PDE、真实参数和合成观测

### 第 81--83 行

```matlab
%% Setup pde and solve
pde = stokes_data_grav_period;
warning('\nChange function to data !!!\n');
```

`stokes_data_grav_period` 返回包含函数句柄的结构体，例如体力、边界数据和精确解。MATLAB 调用无输入函数时，括号可以省略。

`warning` 输出警告。这里的文字是开发阶段遗留提醒，不表示程序发生了运行错误。

### 第 84 行

空行。

### 第 85--88 行

```matlab
xbot = [0:h:1-h h/2:h:1]';
ybot = 0 - slope * xbot;
pt_bot = [xbot, ybot];
```

`xbot` 的排列是

```text
[周期底边顶点坐标; 底边中点坐标]
```

总长度为 $2N_m$。第一个区间不包含 $x=1$，因为它与 $x=0$ 周期等价；第二个区间存放每条边的中点。

底边满足

$$
y=-sx.
$$

`[xbot,ybot]` 水平拼接成二维坐标矩阵。`pt_bot` 后续只用于预先计算边界数据。

### 第 89--93 行

```matlab
xtop = xbot;
ytop = 0.5 - slope * xbot;
pt_top = [xbot, ytop];
```

顶部与底部拥有相同的周期水平坐标，并满足

$$
y=0.5-sx.
$$

### 第 94--98 行

```matlab
% bot
pde.g_N = pde.g_N(pt_top);
pde.g_RN = pde.g_RN(pt_bot);
pde.g_Dn = pde.g_Dn(pt_bot);
```

原本的字段是函数句柄。这里立即在 P2 边界自由度坐标上求值，把它们替换为数值数组。

- `g_N`：顶部 Neumann 数据；
- `g_RN`：底部 Robin 右端；
- `g_Dn`：底部法向速度。

注释 `% bot` 不完全准确，因为第一条赋值处理的是顶部。

### 第 99--104 行

```matlab
pde.exactp = [];
pde.exactux = [];
pde.exactuy = [];
pde.exactu = [];
pde.g_D = [];
```

`[]` 是空数组。清空这些字段表示：

- 不再使用精确解做误差计算；
- 不施加 Dirichlet 速度边界；
- 只保留顶部自然边界、底部滑移边界和周期条件。

### 第 105--108 行

注释说明即将设置真实参数，随后输出 `Variable m0`。

### 第 109 行

```matlab
pde.g_R = 1 + mp('0.1') * ...
    cos(2 * xbot * mp('pi') + 0.1 * pi);
```

构造真实 Robin 系数

$$
m_0(x)
=
1+0.1\cos(2\pi x+0.1\pi).
$$

`cos` 对向量逐元素计算。标量乘向量可直接使用 `*`；若两侧都是同尺寸数组，逐元素乘法应使用 `.*`。

这里混用了多精度 `mp('pi')` 和双精度 `pi`，风格并不统一。

### 第 110--111 行

```matlab
pde.g_R = linearize_bot(pde.g_R);
m0 = pde.g_R;
```

`linearize_bot` 用相邻底部顶点值的平均数覆盖 P2 边中点值：

$$
m_{i+1/2}
=
\frac{m_i+m_{i+1}}{2},
$$

周期末端使用

$$
m_{N_m+1}=m_1.
$$

`m0` 保存真实参数，包括顶点值和边中点值。

### 第 112--113 行

空行。

### 第 114--117 行

```matlab
pde.g_N = linearize_top(pde.g_N);
pde.g_R = linearize_bot(pde.g_R);
pde.g_RN = (pde.g_RN);
pde.g_Dn = (pde.g_Dn);
```

- `linearize_top` 当前实现直接返回输入，所以第一行没有改变数据；
- 第二行重复执行底部线性插值，也是冗余的；
- 圆括号不会改变 `g_RN` 和 `g_Dn`，后两行是无操作语句。

### 第 118--120 行

空行后调用正问题求解器：

```matlab
[soln,eqn,info] = StokesP2P1_periodic(...
    node,elem,bdFlag,pde,option);
```

返回：

- `soln`：速度和压力解；
- `eqn`：离散矩阵、右端和自由度信息；
- `info`：求解器信息。

### 第 121--123 行

```matlab
uh = soln.u;
ph = soln.p;
u_obs = [uh(IUxTop), uh(IUyTop)];
```

把真实参数对应的状态解记为 `uh`、`ph`，并提取顶部两个速度分量作为无噪声合成观测。

后续目标函数只使用 `u_obs(:,1)`，即水平速度；竖直观测虽然保存了，但未参与反演。

## 12. 第 126--161 行：绘制真实解、设置初值并求初始状态

### 第 126 行

```matlab
plt1 = subplot(plot_m, plot_n, 1);
```

选择 $2\times3$ 子图布局中的第一个子图，并保存图形坐标轴句柄。`plt1` 后面没有继续使用。

### 第 127--128 行

```matlab
trisurf(elem,node(:,1),node(:,2),uh(IUxNode),...
        'FaceColor','interp','EdgeColor','interp');
```

在三角网格上绘制顶点水平速度：

- `elem` 给出三角形连接；
- `node(:,1:2)` 给出平面坐标；
- `uh(IUxNode)` 作为高度或颜色数据；
- `'interp'` 表示面颜色和边颜色采用插值。

### 第 129--133 行

```matlab
axis equal;
axis tight;
colorbar;
title('sol u', 'FontSize', 14);
view(2);
```

- `axis equal`：两个坐标轴使用相同比例；
- `axis tight`：坐标范围贴合数据；
- `colorbar`：添加颜色条；
- `title`：设置标题和字体大小；
- `view(2)`：使用二维俯视图。

### 第 134 行

```matlab
stokes_plot_solution;
```

运行另一个脚本，进一步绘制解和边界速度。它依赖当前工作区中的变量，因此与主脚本耦合较强。

### 第 135--137 行

空行后开始设置反演初值。

### 第 138--139 行

```matlab
m = m0 + mp('0.1') * ...
    (sin(xbot * mp('pi') * two) + 0.25);
m = linearize_bot(m);
```

初始猜测为

$$
m^{(0)}
=
m_0+0.1\left[\sin(2\pi x)+0.25\right].
$$

这意味着初值是人为扰动后的真实参数。随后再次保证边中点值等于相邻顶点平均。

### 第 140--146 行

```matlab
pde.g_R = m;
plt3 = subplot(plot_m, plot_n, 3);
plot(...);
legend('m0', 'm');
title('m', 'FontSize', 14);
```

把 PDE 中的 Robin 系数替换为当前猜测，并绘制真实参数和初始参数。

绘图横坐标使用

```matlab
[sft(xbot); sft(xbot)+1]
```

把一个周期复制到 $[1,2]$，从而连续显示两个周期。纵坐标也复制一遍。

`...` 是 MATLAB 续行符，表示下一物理行仍属于同一条语句。

### 第 147--151 行

开始求当前参数下的正问题：

```matlab
[soln,eqn,info] = StokesP2P1_periodic(...);
um = soln.u;
```

数学上得到

$$
\boldsymbol u_m=\boldsymbol u(m).
$$

### 第 152--161 行

绘制当前解和真实解的水平速度差

$$
u_x(m)-u_x(m_0).
$$

第 154 行返回的 `plt2` 句柄后续未使用。其余绘图命令与第 127--133 行相同。

## 13. 第 163--203 行：目标函数导数和第一次伴随梯度

### 第 163--164 行

```matlab
%% Objective
dXidu = two*([um(IUxTop) - u_obs(:,1), xtop*0]);
```

构造目标泛函关于顶部速度的导数：

$$
\Xi_u
=
2
\begin{pmatrix}
u_x-u_x^{\mathrm{obs}}\\
0
\end{pmatrix}.
$$

`xtop*0` 是一种生成同尺寸零向量的写法。

### 第 165--167 行

```matlab
plt4 = subplot(...);
size(xbot);
size(dXidu(:,1));
```

选择第 4 个子图。带分号的 `size(...)` 虽然会计算尺寸，但不显示结果，也不保存结果，因此第 166--167 行没有可观察效果，是调试残留。

### 第 168--171 行

绘制顶部水平速度失配导数，并设置横轴标签、图例和标题。

### 第 172--176 行

空行及注释，说明接下来使用伴随法计算参数梯度。

### 第 177 行

```matlab
pde_adj = pde;
```

MATLAB 结构体采用值语义。这里复制 PDE 数据，之后修改 `pde_adj` 不会修改 `pde`。

### 第 178--179 行

```matlab
pde_adj.f = 0;
pde_adj.fp = 0;
```

伴随问题内部没有体力，也没有非零散度右端。

### 第 180 行

```matlab
pde_adj.g_N = ...
    -[linearize_top(dXidu(:,1)), xtop*0];
```

在顶部施加伴随 Neumann 数据

$$
-\Xi_u
=
-2
\begin{pmatrix}
u_x-u_x^{\mathrm{obs}}\\
0
\end{pmatrix}.
$$

负号来自伴随方程约定

$$
L_u^T u^*=-\Xi_u^T.
$$

### 第 181 行

```matlab
pde_adj.g_R = m;
```

伴随算子使用与正问题相同的底部 Robin 系数。

### 第 182--183 行

```matlab
pde_adj.g_RN = [xbot*0, xbot*0];
pde_adj.g_Dn = [xbot*0, xbot*0];
```

伴随问题的底部 Robin 右端和法向约束数据设为零。

`g_Dn` 在正问题中通常是标量法向速度数据，而这里传入两列零。由于全为零，求解器相关索引仍能取得零值，但数据形状并不清晰，建议重构时统一接口。

### 第 184--185 行

```matlab
[soln,eqn,info] = StokesP2P1_periodic(...);
ustar = soln.u;
```

求解主伴随问题，得到伴随速度 $\boldsymbol u^*$。

### 第 186--190 行

注释说明参数导数来自

$$
\langle u^*,L_m\,\delta m\rangle.
$$

随后

```matlab
dLdm = [um(IUxBot), um(IUyBot)];
```

提取正问题速度在底部的两个分量。

第 190 行以 `%` 开头，因此旧的简单节点乘积公式不会执行。

### 第 191 行

```matlab
dXidm = zeros(Nm, 1);
```

预分配 $N_m\times1$ 梯度向量。预分配可以避免循环中反复扩大数组。

### 第 192--198 行

```matlab
for ii = 1:Nm
    m1 = extend_mid(EI(:, ii));
    dXidm(ii) = integral_robin_P2(...);
end
```

对每个参数基方向：

1. `EI(:,ii)` 取第 $i$ 个参数标准基；
2. `extend_mid` 把周期 P1 基函数扩展到 P2 边界节点；
3. `[m1,m1]` 把同一个标量基函数复制到两个速度分量；
4. `integral_robin_P2` 计算

   $$
   g_i
   =
   \int_{\Gamma_b}
   \psi_i
   \left(
   u_x^*u_x+u_y^*u_y
   \right)\,ds.
   $$

`for ... end` 是 MATLAB 的计数循环。

### 第 199 行

```matlab
dXidm = extend_mid(dXidm);
```

把仅含独立顶点分量的梯度扩展到 P2 边界表示。这一步主要用于绘图；进入反演循环后，梯度仍保持 `Nm × 1`。

### 第 200--203 行

选择第 5 个子图，按几何顺序绘制梯度。`'-bx'` 表示蓝色实线并使用 `x` 标记。

## 14. 第 206--247 行：反演循环、正问题和伴随梯度

### 第 206--209 行

```matlab
%% Iterate to minimize Xi
figure(2);
for k = 1 : max_iteration
```

打开第二个图窗，并执行最多 `max_iteration` 次反演迭代。冒号两侧的空格不影响语义。

### 第 210--214 行

注释说明第一步是求状态。随后：

```matlab
pde_test = pde;
pde_test.g_R = m;
```

复制 PDE 数据并把 Robin 系数替换为当前参数。

### 第 215--217 行

```matlab
[soln,eqn,info] = StokesP2P1_periodic(...);
um = soln.u;
eqn0 = eqn;
```

求当前状态 $\boldsymbol u(m)$。`eqn0` 保存离散方程结构，但当前脚本后面没有使用它。

### 第 218--221 行

```matlab
fprintf(...,
    integral_neumann_P2(...),
    norm(m-m0,Inf),
    norm(um(IUxTop)-u_obs(:,1),Inf));
```

每轮输出：

1. 迭代编号；
2. 目标函数

   $$
   \Xi(m)
   =
   \int_{\Gamma_t}
   (u_x-u_x^{\mathrm{obs}})^2\,ds;
   $$

3. 参数最大范数误差；
4. 顶部水平速度最大范数误差。

`'repeat'` 告诉 `integral_neumann_P2` 把第一个输入同时当作第二个输入，因此计算的是平方积分。

`Inf` 作为 `norm` 的第二参数表示无穷范数：

$$
\|\boldsymbol x\|_\infty=\max_i|x_i|.
$$

### 第 222--224 行

注释说明使用“伴随 + CG”，然后重新分配梯度向量。

实际上调用的是 `cgs`，即 conjugate gradients squared，而不是标准 `cg` 或 `pcg`。

### 第 225--228 行

注释写出主伴随方程

$$
(L_u)^T u^s=-\Xi_u^T.
$$

随后重新计算当前迭代点的

$$
\Xi_u
=
2(u_x-u_x^{\mathrm{obs}},0)^T.
$$

### 第 229--237 行

与第 177--185 行相同：

- 清空伴随体力；
- 把顶部速度失配作为负 Neumann 数据；
- 使用当前参数 `m` 作为底部 Robin 系数；
- 求出伴随速度 `us1`。

这里 `pde_adj=pde` 而不是 `pde_test`，但随后显式设置了与伴随问题有关的字段，因此结果仍使用当前 `m`。

### 第 238--247 行

与第 189--198 行相同，逐个底部参数基函数计算梯度

$$
\left(\nabla\Xi(m)\right)_i
=
\int_{\Gamma_b}
\psi_i\,\boldsymbol u^s\cdot\boldsymbol u\,ds.
$$

这一次没有执行 `extend_mid(dXidm)`，因此 `dXidm` 保持为 $N_m$ 维独立参数梯度。

## 15. 第 249--279 行：稳定矩阵和 Hessian 线性系统

### 第 249--250 行

```matlab
% Stab Matrix
get_robin_stab_mat;
```

`get_robin_stab_mat.m` 是脚本而不是函数，因此它直接在当前工作区创建 `Mstab`。

当前启用分支构造周期一维差分矩阵

$$
M_{\mathrm{stab}}
\approx
\int_{\Gamma_b}
m'(x)v'(x)\,ds.
$$

其矩阵形式包含周期首尾耦合。

### 第 251--252 行

```matlab
gamma_stab = 1e-11;
```

设置很小的稳定化参数 $\gamma=10^{-11}$。

### 第 253--262 行

```matlab
stokes_info.nx = nx;
...
stokes_info.gamma_stab = gamma_stab;
```

把 Hessian-vector product 所需的数据打包到结构体：

- 网格维数；
- 顶部、底部坐标；
- 顶部和底部速度索引；
- 稳定矩阵及系数。

这种方式比继续增加函数位置参数更清晰。

### 第 263--264 行

```matlab
dXidm_stab = dXidm;
```

变量名暗示它应包含正则化梯度，但当前代码只是复制数据梯度。

如果目标函数真的包含

$$
\frac{\gamma}{2}m^TM_{\mathrm{stab}}m,
$$

则梯度应为

$$
\nabla\Xi_{\mathrm{reg}}
=
\nabla\Xi
+
\gamma M_{\mathrm{stab}}m.
$$

当前主程序没有加上第二项，而 `stokes_hessian.m` 却给 Hessian-vector product 加了

$$
\gamma M_{\mathrm{stab}}\delta m.
$$

因此这里更准确地说是“只在更新方程左端加入平滑阻尼”，并不是一个完全一致的 Tikhonov 正则化 Newton 系统。

### 第 265--266 行

空行后重新构造单位矩阵。`EI` 之前已经存在，所以这一步是冗余的。

### 第 267 行

```matlab
if true
```

条件恒为真，因此程序始终进入迭代矩阵自由求解分支。第 272--279 行的显式矩阵分支不会执行，除非手动把这里改为 `false`。

### 第 268--271 行

```matlab
[dm, flg, relres, niter, resvec] = ...
    cgs(@(m1) stokes_hessian(...,m1,...), ...
        dXidm_stab(:), 1e-10, 50);
```

`cgs` 求解

$$
H\,dm=g.
$$

第一个输入不是矩阵，而是匿名函数：

```matlab
@(m1) stokes_hessian(...,m1,...)
```

它接收任意方向 `m1`，返回 $Hm_1$，因此不需要显式存储 Hessian。

其他参数为：

- 右端 `dXidm_stab(:)`；`(:)` 强制转成列向量；
- 相对残差容限 $10^{-10}$；
- 最大 50 次 Krylov 迭代。

输出为：

- `dm`：参数更新方向；
- `flg`：收敛标志；
- `relres`：最终相对残差；
- `niter`：迭代次数；
- `resvec`：残差历史。

脚本只打印 `niter` 和 `relres`，没有检查 `flg`。稳健实现中应检查 `flg` 是否为零。

### 第 272--279 行

这是不会执行的备选分支：

```matlab
d2Xidm2 = zeros(Nm);
for i = 1:Nm
    d2Xidm2(:,i) = stokes_hessian(...,EI(:,i),...);
end
dm = d2Xidm2 \ dXidm_stab(:);
```

它逐列作用 Hessian 来显式组装矩阵：

$$
H
=
\begin{bmatrix}
He_1&He_2&\cdots&He_{N_m}
\end{bmatrix},
$$

然后使用反斜杠求解 $Hdm=g$。

MATLAB 的 `\` 表示线性系统求解，不应写成 `inv(H)*g`。

## 16. 第 281--294 行：显示更新并修改参数

### 第 281 行

```matlab
if k <= 4
```

第二个图窗只安排了四个纵向子图，因此只绘制前四次更新。

### 第 282 行

```matlab
subplot(4,1,k)
```

选择 $4\times1$ 布局中的第 `k` 个子图。末尾没有分号，某些 MATLAB 版本可能显示返回的图形句柄。

### 第 283--284 行

绘制：

- 蓝色圆点：求得的独立参数更新 `dm`；
- 红色叉号：当前参数误差 $m-m_0$。

如果 Hessian 模型准确，则希望

$$
dm\approx m-m_0,
$$

因为更新使用 $m\leftarrow m-dm$。

`dm` 只有 $N_m$ 个顶点值，所以直接复制两个周期；`m-m0` 已经是 P2 边界布局，因此先重排再复制。

### 第 285--287 行

设置图例，并使用

```matlab
s = sprintf('iter :%d', k);
```

生成包含迭代编号的标题字符串。`sprintf` 返回字符串但不直接打印。

### 第 288 行

结束绘图 `if` 块。

### 第 289--292 行

```matlab
% update m
mbefore = m;
m = mbefore - extend_mid(dm);
```

保存更新前参数，并执行

$$
m^{(k+1)}
=
m^{(k)}-\delta m.
$$

`extend_mid(dm)` 把独立周期顶点更新扩展到底部 P2 节点和边中点。

代码没有线搜索、信赖域、参数正值约束或停止判据，因此固定执行五次完整 Newton/Gauss--Newton 步。如果步长过大，可能导致 $m$ 变成负值或目标函数增大。

### 第 293 行

空行。

### 第 294 行

```matlab
end
```

结束反演循环。

### 第 295--296 行

文件末尾空行。

## 17. `stokes_hessian.m` 实际做了什么

主脚本最关键的外部函数是 `stokes_hessian.m`。它不是返回 Hessian 矩阵，而是返回 Hessian 对输入方向 `m1` 的作用。

### 17.1 参数方向的状态增量

代码设置

```matlab
pde_du.g_R = linearize_bot(m);
pde_du.g_RN = [-um(IUxBot), m1, ...
               -um(IUyBot), m1];
```

`StokesP2P1_periodic` 对四列 `g_RN` 的特殊解释是先逐列插值，再计算

```matlab
gp(:,1) = gp(:,1).*gp(:,2);
gp(:,2) = gp(:,3).*gp(:,4);
```

因此底部右端是

$$
-m_1
\begin{pmatrix}
u_x\\
u_y
\end{pmatrix}.
$$

这对应线性化 Robin 项

$$
L_u\,\delta\boldsymbol u
=
-L_m m_1.
$$

### 17.2 增量伴随问题

状态增量 `du1` 在顶部产生观测增量。函数随后求解

```matlab
pde_adj3.g_N = ...
    -[2*linearize_top(du1(IUxTop)), xtop*0];
```

即

$$
L_u^T\boldsymbol u_3^*
=
-2C^TC\,\delta\boldsymbol u.
$$

### 17.3 投影回参数空间

对每个参数基函数 $\psi_i$，计算

$$
\int_{\Gamma_b}
\psi_i\,
\boldsymbol u_3^*\cdot\boldsymbol u\,ds.
$$

最后加上

$$
\gamma M_{\mathrm{stab}}m_1.
$$

因此该函数实现的是 Gauss--Newton 型 Hessian-vector product，而不是包含所有二阶状态导数项的完整 Newton Hessian。

## 18. 主要辅助函数

### `integral_robin_P2`

在所有 Robin 边上使用六阶一维求积，计算

$$
\int_{\Gamma_b}
f\,g\,ds
$$

或

$$
\int_{\Gamma_b}
f\,g\,h\,ds.
$$

如果输入有两列，则先对两个速度分量求和。

### `integral_neumann_P2`

与前者结构相同，但积分区域是标记为 Neumann 的顶部边界。主脚本用它计算

$$
\int_{\Gamma_t}
(u_x-u_x^{\mathrm{obs}})^2\,ds.
$$

### `extend_mid`

输入 $N_m$ 个周期顶点值，输出 $2N_m$ 个 P2 边界值：

$$
\operatorname{extend\_mid}(m)
=
\begin{bmatrix}
m_1,\ldots,m_{N_m},
\frac{m_1+m_2}{2},\ldots,
\frac{m_{N_m}+m_1}{2}
\end{bmatrix}^T.
$$

### `linearize_bot`

输入已经具有 $2N_m$ 个分量的数组，保留前半部分顶点值，用相邻顶点平均覆盖后半部分边中点值。

### `linearize_top`

当前启用分支只是返回输入。函数中保留了一个被 `if 0` 禁用的顶点到边中点插值实现。

### `get_robin_stab_mat`

创建周期一维刚度型矩阵。当前分支等价于周期有限差分 Laplace 算子的缩放版本。

## 19. 算法流程总结

整个脚本可以概括为：

1. 在参考矩形上生成网格并标记边界；
2. 将网格倾斜成平行四边形；
3. 构造空间变化的真实底部 Robin 系数 $m_0$；
4. 求解线性 Stokes 正问题，生成顶部速度观测；
5. 构造扰动后的初始参数 $m^{(0)}$；
6. 对每次反演迭代：
   1. 求当前正问题；
   2. 计算顶部水平速度失配；
   3. 求主伴随问题；
   4. 计算参数梯度；
   5. 通过 Hessian-vector product 和 `cgs` 求 Gauss--Newton 步；
   6. 更新底部参数。

用伪代码表示为

```text
生成观测 u_obs = C u(m0)
选择初值 m

for k = 1,...,5
    解状态方程       L(m) u = f
    计算残差         r = C u - u_obs
    解伴随方程       L(m)^T u* = -2 C^T r
    计算梯度         g_i = ∫Γb ψ_i u*·u ds
    用 CGS 解         H_GN dm = g
    更新参数         m = m - dm
end
```

## 20. 阅读和使用这份代码时应注意的问题

### 20.1 这是线性 Stokes 反演

状态方程的黏度不依赖速度，应力--应变关系是线性的。它没有调用 `NonlinearStokesP2P1_periodic.m`，因此不是 Glen 非线性黏度反演。

### 20.2 变量 `m` 是底部 Robin 系数

它更适合命名为 `beta`。如果把它与 Weertman 指数或 Glen 指数也记作 $m$，很容易混淆。

### 20.3 只使用水平速度观测

`u_obs` 保存两个速度分量，但目标泛函的竖直分量被显式设为零。若希望同时使用两个分量，需要修改 `dXidu`、目标积分和 Hessian 中的顶部右端。

### 20.4 观测无噪声且来自同一个离散模型

这属于 inverse crime：生成观测和执行反演使用完全相同的网格、有限元空间和求解器。因此结果通常比真实数据反演乐观。

更严格的数值实验应使用：

- 更细网格生成观测；
- 不同网格执行反演；
- 加入可控噪声；
- 根据噪声水平选择正则化参数。

### 20.5 正则化并不完全一致

Hessian-vector product 包含

$$
\gamma M_{\mathrm{stab}}\delta m,
$$

但梯度没有包含

$$
\gamma M_{\mathrm{stab}}m.
$$

因此代码并没有严格求解某个包含标准二次正则项的目标泛函。

### 20.6 没有全局化策略

更新没有：

- 线搜索；
- Levenberg--Marquardt 阻尼；
- 信赖域；
- 参数上下界；
- 正值参数化；
- 基于梯度或目标函数的停止条件。

这意味着它更适合作为算法推导和小规模合成实验，而不是直接用于真实冰盖数据。

### 20.7 `cgs` 不等于标准共轭梯度

若 Gauss--Newton Hessian 确实对称正定，通常更自然的选择是 `pcg`。当前代码使用 `cgs`，它可以处理更一般的线性系统，但应检查返回标志 `flg`。

## 21. 一句话理解这份程序

`codesinmypaper.m` 用 P2--P1 有限元求解倾斜周期区域上的线性 Stokes 方程，把顶部水平速度误差通过伴随方程传回底部，再用矩阵自由 Gauss--Newton 方法更新空间变化的底部 Robin 摩擦系数。
