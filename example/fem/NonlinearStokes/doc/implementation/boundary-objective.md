# `NonlinearStokesAdjointInversionBoundaryObjective.m` 详解

本文解释脚本

```text
example/fem/NonlinearStokes/NonlinearStokesAdjointInversionBoundaryObjective.m
```

的目的、数学含义、代码结构和主要输出。

这个脚本做的是一个二维非线性 Stokes 冰流模型中的底部滑移参数反演问题：

- 正问题：给定底部滑移系数 $\beta(x)$，求速度 $u$ 和压力 $p$；
- 观测：取顶部边界上的水平速度；
- 反问题：根据顶部速度观测恢复底部滑移系数 $\beta(x)$；
- 优化变量：不是直接用 $\beta$，而是用

$$
q=\log(\beta).
$$

使用 $q$ 的直接好处是保证

$$
\beta=\exp(q)>0.
$$

脚本的核心算法是：

1. 用阻尼 Picard 迭代求非线性 Stokes 正问题；
2. 用一致切线矩阵解伴随方程，得到目标函数梯度；
3. 用增量状态方程和增量伴随方程构造矩阵自由的 Gauss--Newton Hessian-vector product；
4. 用 PCG 解阻尼 Gauss--Newton / Levenberg--Marquardt 步；
5. 用回溯线搜索接受下降步。

## 1. 脚本整体目标

脚本开头注释说明它是从 `NSAdjointInversion.m` 修改而来，主要区别是目标函数换成了顶部边界积分误差：

$$
J(q)
=
\frac{1}{2}
\frac{
\int_{\Gamma_t} (u(q)-u_{\mathrm{obs}})^2\,ds
}{
\int_{\Gamma_t} u_{\mathrm{obs}}^2\,ds
}.
$$

这里：

- $\Gamma_t$ 是顶部边界；
- $u(q)$ 是给定 $q$ 后求出的正问题速度；
- $u_{\mathrm{obs}}$ 是合成观测数据；
- 分母用于归一化，使目标函数变成无量纲量；
- 脚本中没有额外正则项。

代码中对应的是：

```matlab
residual = u(topDof)-dataObs;
dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
objective = dataObjective;
```

其中：

```matlab
dataNormSquared = max(topWeight'*(dataObs.^2),eps);
```

表示离散后的

$$
\int_{\Gamma_t} u_{\mathrm{obs}}^2\,ds.
$$

## 2. 几何、网格和边界条件

脚本使用一个二维矩形区域作为初始网格：

```matlab
L = 1;
H = 0.5;
slope = 0.1;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
```

然后设置边界：

```matlab
topBoundaryExpression = sprintf('y==%.17g',H);
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
```

含义是：

- 顶部 `y==H` 是 Neumann 边界；
- 底部 `y==0` 是 Robin 边界，也就是滑移边界；
- 后续通过坐标变换把矩形变成倾斜 slab：

```matlab
node(:,2) = node(:,2)-slope*node(:,1);
```

因此实际顶部边界为

$$
y=H-\mathrm{slope}\cdot x.
$$

脚本还启用了周期边界：

```matlab
option.periodic = true;
option.periodic_x = [0,L];
```

所以 $x=0$ 和 $x=L$ 两侧被识别为周期边界。

## 3. 顶部观测自由度

速度使用 P2 元，因此速度自由度不仅包括网格顶点，也包括边中点。

代码先构造所有 P2 速度自由度坐标：

```matlab
[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
```

然后找出顶部边界上的速度自由度：

```matlab
surfaceLevel = H-slope*uNode(:,1);
topDof = find(abs(uNode(:,2)-surfaceLevel)<tolGeometry ...
            & uNode(:,1)<L-tolGeometry);
```

这里有一个条件：

```matlab
uNode(:,1)<L-tolGeometry
```

它排除了周期右端点，避免 $x=0$ 和 $x=L$ 两侧重复计数。

随后按 $x$ 坐标排序：

```matlab
[~,order] = sort(uNode(topDof,1));
topDof = topDof(order);
xObs = uNode(topDof,1);
```

这些 `topDof` 就是目标函数中使用的顶部速度观测自由度。

## 4. 顶部边界积分权重

目标函数是边界积分形式。脚本用 `boundaryweights` 构造离散权重：

```matlab
topWeight = boundaryweights(xObs,L,slope);
```

局部函数为：

```matlab
function weight = boundaryweights(xObs,L,slope)
    nObs = numel(xObs);
    assert(nObs > 0,'No top-boundary observation dofs were found.');
    weight = sqrt(1+slope^2)*(L/nObs)*ones(nObs,1);
end
```

由于顶部边界是倾斜直线

$$
y=H-\mathrm{slope}\cdot x,
$$

其弧长微元为

$$
ds=\sqrt{1+\mathrm{slope}^2}\,dx.
$$

所以每个观测点分配的权重近似为

$$
w_i=\sqrt{1+\mathrm{slope}^2}\frac{L}{N_{\mathrm{obs}}}.
$$

于是

```matlab
topWeight'*(residual.^2)
```

近似

$$
\int_{\Gamma_t} (u-u_{\mathrm{obs}})^2\,ds.
$$

## 5. 非线性 Stokes 模型设置

脚本中的 PDE 参数为：

```matlab
pde.A = 1;
pde.n = 3;
pde.m = 1/3;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];
```

其中：

- `pde.A` 是 Glen 定律中的流动率；
- `pde.n = 3` 是 Glen 指数；
- `pde.m = 1/3` 是 Weertman 滑移律指数；
- `pde.rho` 和 `pde.gravity` 给出重力体力；
- `pde.g_N = []` 表示顶部没有额外给定牵引。

正问题由

```matlab
NonlinearStokesP2P1(node,elem,bdFlag,pde,option)
```

求解。该函数使用：

- P2 速度；
- P1 压力；
- 阻尼 Picard 迭代处理非线性黏度和非线性滑移；
- 可选的一致切线矩阵 `eqn.tangent`。

脚本中的正问题选项为：

```matlab
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
option.assemble_tangent = true;
```

其中最重要的是：

- `option.tol = 1e-11`：Picard 相对变化停止阈值；
- `option.damping = 0.8`：Picard 阻尼系数；
- `option.assemble_tangent = true`：求解完成后装配一致切线矩阵，用于伴随和增量方程。

## 6. 参数化方式

底部滑移系数 $\beta(x)$ 使用周期 P1 参数化。

代码为：

```matlab
Nm = round(L/h);
xBeta = (0:Nm-1)'*h;
betaTrue = 1+0.1*cos(2*pi*xBeta/L);
betaInitial = betaTrue+0.1*(sin(2*pi*xBeta/L)+0.25);
qTrue = log(betaTrue);
q = log(betaInitial);
```

当前设置下：

```text
L = 1
h = 1/8
Nm = 8
```

所以有 8 个参数自由度。

`periodicP1` 用于把这些节点值插值成周期函数：

```matlab
function value = periodicP1(x,xNode,nodalValue,L)
    xWrapped = mod(x,L);
    value = interp1([xNode;L],[nodalValue;nodalValue(1)],...
        xWrapped,'linear');
end
```

这表示：

- 参数在 $[0,L)$ 上有 `Nm` 个节点；
- $x=L$ 处的值等于 $x=0$ 处的值；
- 插值是线性的；
- 超出周期区间时用 `mod(x,L)` 包回。

在正问题中，脚本通过

```matlab
beta = exp(q(:));
pde.beta = @(pt) periodicP1(pt(:,1),xBeta,beta,L);
```

把参数向量 `q` 转换成底部滑移函数。

## 7. 合成观测数据

脚本不是读取外部观测数据，而是用真实参数 `qTrue` 生成合成观测：

```matlab
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
```

也就是说：

1. 用 `qTrue` 解一次正问题；
2. 取顶部自由度上的水平速度；
3. 把它作为观测数据 `dataObs`。

这种方式适合做算法验证，因为真实参数 `betaTrue` 已知，可以检查恢复误差。

归一化因子为：

```matlab
dataNormSquared = max(topWeight'*(dataObs.^2),eps);
```

如果观测速度非常小，`max(...,eps)` 防止除以零。

## 8. 反演历史量

脚本预分配了一组历史量：

```matlab
history.objective
history.dataResidual
history.parameterError
history.parameterErrorLinf
history.parameterErrorRelativeLinf
history.gradientNorm
history.picardSteps
```

含义分别是：

- `objective`：目标函数值；
- `dataResidual`：归一化数据残差范数；
- `parameterError`：$\beta$ 的相对 $L^2$ 向量误差；
- `parameterErrorLinf`：$\beta$ 的绝对 $\ell^\infty$ 误差；
- `parameterErrorRelativeLinf`：相对 $\ell^\infty$ 误差；
- `gradientNorm`：梯度向量二范数；
- `picardSteps`：本轮正问题 Picard 迭代步数。

这些量只用于诊断和画图，不参与算法本身。

## 9. 每轮反演迭代

主循环为：

```matlab
for k = 1:maxInverseIt
    ...
end
```

当前最大迭代数为：

```matlab
maxInverseIt = 10;
```

每一轮包含以下步骤。

### 9.1 解正问题

```matlab
[u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
    node,elem,bdFlag,xBeta,L);
```

这里：

- `q` 是当前参数；
- `uWarm` 是上一轮速度解，用作初值；
- `u` 是当前正问题速度；
- `eqn` 包含离散矩阵、切线矩阵和参数导数接口；
- `forwardInfo.itStep` 是 Picard 步数；
- `forwardInfo.converged` 表示正问题是否收敛。

脚本要求正问题必须收敛：

```matlab
assert(forwardInfo.converged,...
    'Forward solve failed at inverse iteration %d.',k);
```

### 9.2 计算残差和目标函数

```matlab
residual = u(topDof)-dataObs;
dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
objective = dataObjective;
```

这就是离散化后的

$$
J(q)
=
\frac{1}{2}
\frac{
\sum_i w_i (u_i-u_{\mathrm{obs},i})^2
}{
\sum_i w_i u_{\mathrm{obs},i}^2
}.
$$

因为没有正则项，所以

```matlab
objective = dataObjective;
```

### 9.3 装配参数导数矩阵 `G`

```matlab
G = assembleparameterderivative(eqn,q,xBeta,L,Nm);
```

`G` 的含义是离散残差对参数 `q` 的导数：

$$
G = \frac{\partial R}{\partial q}.
$$

其中正问题残差写成

$$
R(U,q)=0.
$$

这里

$$
U=
\begin{bmatrix}
u\\
p\\
\lambda
\end{bmatrix}
$$

包含速度、压力和约束乘子。

局部函数中逐列构造 `G`：

```matlab
for j = 1:Nm
    direction = zeros(Nm,1);
    direction(j) = 1;
    deltaBeta = beta.*direction;
    directionFunction = @(pt) periodicP1(...
        pt(:,1),xBeta,deltaBeta,L);
    G(:,j) = eqn.applyBetaDerivative(directionFunction);
end
```

注意这里有链式法则：

$$
\beta=\exp(q),
\qquad
\delta\beta=\beta\,\delta q.
$$

所以第 `j` 个方向上，

```matlab
deltaBeta = beta.*direction;
```

表示对 $q_j$ 求导时对应的 $\delta\beta$。

### 9.4 构造目标函数对状态变量的导数

```matlab
observationGradient = zeros(size(eqn.tangent,1),1);
observationGradient(topDof) = topWeight.*residual/dataNormSquared;
```

目标函数对观测速度的导数为

$$
\frac{\partial J}{\partial u_i}
=
\frac{w_i(u_i-u_{\mathrm{obs},i})}
{\sum_j w_j u_{\mathrm{obs},j}^2}.
$$

这正是：

```matlab
topWeight.*residual/dataNormSquared
```

由于目标函数只依赖顶部水平速度，所以只有 `topDof` 上的分量非零。

### 9.5 解伴随方程

```matlab
adjoint = eqn.tangent'\(-observationGradient);
gradient = G'*adjoint;
```

设

$$
A = R_U(U,q)=\frac{\partial R}{\partial U}.
$$

代码中的

```matlab
eqn.tangent
```

就是一致切线矩阵 $A$。

伴随方程为

$$
A^T \lambda = -J_U^T.
$$

代码中：

```matlab
adjoint = eqn.tangent'\(-observationGradient);
```

就是在解这个线性系统。

然后 reduced gradient 为

$$
\nabla_q J
=
G^T\lambda.
$$

代码中：

```matlab
gradient = G'*adjoint;
```

这里没有额外正则项，所以梯度只有这一项。

## 10. 什么是一致切线矩阵

详细概念见 [tangent-matrix.md](tangent-matrix.md)。本脚本只依赖下面这个关系：把收敛后的正问题写成离散残差

$$
R(U,q)=0,
$$

其中 `U` 包含速度、压力和约束乘子，`q=\log\beta` 是反演参数。对状态变量的一阶导数

$$
R_U(U,q)=\frac{\partial R}{\partial U}(U,q).
$$

在代码中就是

```matlab
eqn.tangent
```

### 10.1 为什么反演必须用一致切线矩阵

反演中需要状态对参数的导数。对

$$
R(U(q),q)=0
$$

求导，得到

$$
R_U\,\delta U + R_q\,\delta q = 0.
$$

因此增量状态方程是

$$
R_U\,\delta U = -R_q\,\delta q.
$$

如果这里的 $R_U$ 不是完整残差的 Jacobian，而是 Picard 矩阵，那么求出来的 $\delta U$ 就不是正问题解对参数的正确一阶变化。

伴随梯度也依赖同一个矩阵：

$$
R_U^T z=-J_U^T.
$$

所以如果伴随方程里用错了矩阵，得到的梯度通常会和有限差分梯度不一致。

这就是为什么脚本设置：

```matlab
option.assemble_tangent = true;
```

并在伴随和增量状态方程中使用：

```matlab
eqn.tangent
```

而不是使用 Picard 迭代中的冻结系数矩阵。

### 10.2 在代码中在哪里装配

`eqn.tangent` 由 `NonlinearStokesP2P1.m` 在正问题收敛后装配。

相关代码逻辑是：

```matlab
if option.assemble_tangent
    Kt = assembleviscoustangent(u);
    Kbt = assemblebedtangent(u);
    tangentM = [Kt+Kbt, B'; B, sparse(Np,Np)];
    eqn.tangent = [tangentM, C'; C, sparse(nConstraint,nConstraint)];
    eqn.applyBetaDerivative = @assemblebetadirection;
end
```

其中：

- `assembleviscoustangent(u)` 装配黏性项的一致切线；
- `assemblebedtangent(u)` 装配底部滑移项的一致切线；
- `B` 是速度-压力耦合块；
- `C` 是周期条件、法向约束、压力约束等约束矩阵；
- `eqn.applyBetaDerivative` 用于计算残差对参数方向的导数。

因此 `eqn.tangent` 是包含速度、压力和约束的完整增广 Jacobian。

## 11. 为什么伴随问题没有 Picard 迭代

正问题是非线性的：

$$
R(U,q)=0.
$$

所以需要 Picard 迭代。

但伴随问题是在已经收敛的正问题解 $U$ 上，对线性化矩阵求解：

$$
R_U(U,q)^T\lambda=-J_U^T.
$$

这是一个线性系统。

因此脚本中的伴随求解：

```matlab
eqn.tangent'\(...)
```

不会再进行 Picard 迭代。它只解一次线性方程。

## 12. 和文献中“再解一个伴随 Stokes 问题”的关系

一些冰盖反演或 Stokes 反演文献会把伴随方程写成另一个 Stokes 型 PDE。形式上，它看起来像正问题：

- 同样有伴随速度和伴随压力；
- 同样满足类似 Stokes 的弱形式；
- 但源项、边界条件和观测项换成了由目标函数导数给出的项。

因此文献中经常会说：伴随方程可以通过“修改正问题的源项和边界条件”来求解。

这句话在连续 PDE 层面是合理的。但当前代码选择的是另一种等价的离散实现方式：不把伴随方程重新包装成一个新的 PDE 求解器，而是直接解离散伴随线性系统。

正问题离散残差写成

$$
R(U,q)=0,
$$

其中

$$
U=
\begin{bmatrix}
u\\
p\\
\lambda
\end{bmatrix}.
$$

在当前正问题解处线性化，得到一致切线矩阵

$$
A=R_U(U,q).
$$

伴随方程离散后就是

$$
A^T z=-J_U^T.
$$

代码中：

```matlab
adjoint = eqn.tangent'\(-observationGradient);
```

正是在解这个系统。

也就是说，文献中的说法和代码中的实现可以对应为：

```text
连续 PDE 层面:
  解一个伴随 Stokes 方程，源项和边界条件来自目标函数导数。

离散矩阵层面:
  解一致切线矩阵转置系统 (R_U)' * adjoint = -J_U。
```

这两者目标相同，区别只是实现层次不同。

### 12.1 为什么没有直接复用 `NonlinearStokesP2P1`

一个容易误解的地方是：既然伴随方程像 Stokes 方程，为什么不直接调用一次正问题求解器 `NonlinearStokesP2P1`，只改源项和边界条件？

原因是当前非线性模型中，伴随方程需要的矩阵不是 Picard 迭代中的冻结系数矩阵，而是完整非线性残差的一致切线矩阵的转置：

$$
R_U(U,q)^T.
$$

`NonlinearStokesP2P1` 求正问题时，Picard 迭代每一步解的是冻结黏度、冻结滑移系数后的 Stokes 型线性系统。这个 Picard 矩阵适合稳健地求非线性正问题，但它通常不是完整残差的 Jacobian。

对于 Glen 非线性黏度和非线性 Weertman 滑移律，完整一致切线矩阵还包含：

- 黏度 $\eta(u)$ 对速度 $u$ 的导数；
- 底部滑移牵引对速度 $u$ 的导数；
- 压力块和不可压约束块；
- 周期边界和法向约束对应的离散约束块。

这些项对伴随梯度的一致性很重要。如果用 Picard 矩阵替代 `eqn.tangent`，得到的伴随变量一般就不对应原始非线性残差的精确线性化，梯度也可能和有限差分不一致。

因此脚本先让正问题求解器在收敛解处额外装配一致切线矩阵：

```matlab
option.assemble_tangent = true;
```

然后用

```matlab
eqn.tangent'
```

直接作为伴随矩阵。

### 12.2 如果要写一个显式 `solveadjoint` 函数

也可以把当前代码包装成一个显式的伴随求解器，例如概念上写成：

```matlab
function adjoint = solveadjoint(eqn,observationGradient)
    adjoint = eqn.tangent'\(-observationGradient);
end
```

这样接口上更接近“求解一个伴随问题”，但本质上仍然是在解同一个离散转置系统。

如果进一步想把它写成和正问题类似的 PDE 装配形式，就必须保证装配出来的矩阵与 `eqn.tangent'` 完全一致，包括所有黏度导数、滑移导数和约束块。否则它只是在解一个近似伴随问题。

### 12.3 简短结论

所以当前实现不是“不需要伴随方程”，而是：

```text
当前实现已经在离散层面解了伴随方程。
```

它没有再次调用正问题求解器，是因为：

1. 伴随问题在当前正问题解处是线性的；
2. 所需矩阵已经由 `eqn.tangent` 给出；
3. 直接解 `eqn.tangent' \ rhs` 比重新包装成一个正问题求解器更直接；
4. 这样可以避免误用 Picard 矩阵代替一致切线矩阵。

## 13. Gauss--Newton Hessian-vector product

反演中没有显式装配完整 Hessian，而是用函数句柄表示 Hessian-vector product：

```matlab
hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
    topDof,topWeight,dataNormSquared,lambda);
```

然后用 PCG 解：

```matlab
[step,flag,relativeResidual,pcgIt] = pcg(...
    hessian,-gradient,pcgTolerance,pcgMaxIt);
```

`gaussnewtonproduct` 做的是：

```matlab
incrementalState = eqn.tangent\(-G*direction);
```

这对应增量状态方程：

$$
A\,\delta U = -G\,\delta q.
$$

其中：

- $A=R_U$；
- $G=R_q$；
- `direction` 是参数方向 $\delta q$；
- `incrementalState` 是该参数方向引起的状态变化 $\delta U$。

然后构造增量观测梯度：

```matlab
incrementalObservation(topDof) = ...
    topWeight.*incrementalState(topDof)/dataNormSquared;
```

这相当于把 $\delta u$ 代入目标函数的二阶观测项。

再解增量伴随：

```matlab
incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
```

最后得到 Hessian-vector product：

```matlab
product = G'*incrementalAdjoint+lambda*direction;
```

这对应

$$
\left(H_{\mathrm{GN}}+\lambda I\right)\delta q.
$$

其中：

- $H_{\mathrm{GN}}$ 是 Gauss--Newton Hessian 近似；
- $\lambda I$ 是 Levenberg--Marquardt 阻尼；
- 这里不是完整二阶 Hessian，因为没有包含 $R_{UU}$、$R_{Uq}$、$R_{qq}$ 等二阶项。

## 14. LM 阻尼参数 `lambda`

脚本初始设置：

```matlab
lambda = 1e-7;
```

它用于解阻尼 Gauss--Newton 系统：

$$
\left(H_{\mathrm{GN}}+\lambda I\right)s=-g.
$$

其中：

- $s$ 是参数更新步；
- $g$ 是梯度；
- $\lambda$ 越大，步子越保守；
- $\lambda$ 越小，越接近普通 Gauss--Newton。

如果 line search 成功接受下降步：

```matlab
lambda = max(lambda/3,1e-12);
```

阻尼减小。

如果找不到下降步：

```matlab
lambda = 10*lambda;
```

阻尼增大。

这是一种简单的 Levenberg--Marquardt 调节策略。

## 15. 回溯线搜索

PCG 得到方向 `step` 后，脚本不直接接受，而是做回溯线搜索：

```matlab
accepted = false;
stepLength = 1;
for lineSearchIt = 1:10
    qTrial = q+stepLength*step;
    ...
    if trialObjective < objective
        q = qTrial;
        uWarm = uTrial;
        lambda = max(lambda/3,1e-12);
        accepted = true;
        break
    end
    stepLength = stepLength/2;
end
```

流程是：

1. 先尝试完整步长 `stepLength=1`；
2. 如果目标函数下降，则接受；
3. 如果不下降，则步长减半；
4. 最多尝试 10 次；
5. 如果仍然失败，则不更新 `q`，只增大 `lambda`。

线搜索中每一次 trial 都会调用一次正问题求解器。

## 16. 停止条件

脚本有两个主要停止条件。

### 16.1 梯度足够小

```matlab
if norm(gradient) <= gradientTolerance
    ...
    break
end
```

其中：

```matlab
gradientTolerance = 1e-9;
```

### 16.2 步长足够小

```matlab
if norm(step) <= stepTolerance*max(1,norm(q))
    ...
    break
end
```

其中：

```matlab
stepTolerance = 1e-7;
```

它表示相对参数尺度而言，更新步已经很小。

如果达到最大反演迭代数：

```matlab
maxInverseIt = 10;
```

循环也会结束。

## 17. 输出表格解释

脚本使用 `printiterationheader` 和 `printiterationrow` 输出表格。

表头为：

```text
it objective betaL2rel betaLinfAbs betaLinfRel |grad| fPicard pcgIt pcgRel ls stop
```

各列含义如下。

### `it`

反演迭代编号。

### `objective`

当前目标函数：

$$
\frac{1}{2}
\frac{
\int_{\Gamma_t}(u-u_{\mathrm{obs}})^2\,ds
}{
\int_{\Gamma_t}u_{\mathrm{obs}}^2\,ds
}.
$$

### `betaL2rel`

恢复参数与真实参数之间的相对二范数误差：

```matlab
norm(betaCurrent-betaTrue)/norm(betaTrue)
```

### `betaLinfAbs`

绝对无穷范数误差：

```matlab
norm(betaCurrent-betaTrue,inf)
```

### `betaLinfRel`

相对无穷范数误差：

```matlab
norm(betaCurrent-betaTrue,inf)/norm(betaTrue,inf)
```

### `|grad|`

目标函数对参数 `q` 的梯度二范数：

```matlab
norm(gradient)
```

### `fPicard`

本轮正问题 Picard 迭代步数：

```matlab
forwardInfo.itStep
```

注意它只统计主正问题求解，不统计线搜索 trial 中额外正问题的 Picard 步数。

### `pcgIt`

PCG 求解 Gauss--Newton 步所用迭代次数。

如果因为梯度已经满足停止条件而没有进入 PCG，则显示 `-`。

### `pcgRel`

PCG 返回的相对残差。

### `ls`

本轮 line search 中尝试了多少个 trial 正问题。

如果因为梯度停止而没有 line search，则为 `0`。

### `stop`

停止原因：

- `grad`：梯度范数达到阈值；
- `step`：步长达到阈值；
- `-`：本轮没有停止，继续迭代。

## 18. Summary 输出

循环结束后，脚本输出：

```matlab
fprintf('  optimization forward solves: %d\n',...
    optimizationForwardSolves);
```

`optimizationForwardSolves` 统计优化过程中调用 `solveforward` 的次数，包括：

- 每轮主正问题；
- 第一次迭代导数检查中的两次有限差分正问题；
- line search 中的 trial 正问题。

不包括：

- 生成合成观测数据的 truth solve；
- 最后画图前重新求 recovered 解的后处理 solve。

随后输出最终误差：

```matlab
fprintf('  final beta Linf error: %.04e, relative Linf: %.04e\n',...
    betaErrorLinf,betaErrorRelativeLinf);
```

以及第一次迭代的导数检查：

```matlab
fprintf(['  derivative check: state %.04e, grad %.04e, GN %.04e ',...
         '(FD %.04e, adj %.04e)\n'],...
    derivativeCheck.stateError,derivativeCheck.gradientError,...
    derivativeCheck.gaussNewtonError,derivativeCheck.finiteDifference,...
    derivativeCheck.adjointDirection);
```

## 19. 导数检查 `verifyderivatives`

第一次反演迭代时，脚本调用：

```matlab
derivativeCheck = verifyderivatives(...);
```

它选取一个固定方向：

```matlab
direction = sin((1:numel(q))');
direction = direction/norm(direction);
epsilon = 1e-3;
```

然后比较三个量。

### 19.1 状态增量检查

用中心差分：

$$
\frac{u(q+\epsilon d)-u(q-\epsilon d)}{2\epsilon}
$$

与线性化状态方程给出的

$$
\delta U
$$

比较。

代码为：

```matlab
incrementalState = eqn.tangent\(-G*direction);
finiteDifferenceState = (uPlus-uMinus)/(2*epsilon);
```

### 19.2 梯度方向导数检查

用目标函数中心差分：

$$
\frac{J(q+\epsilon d)-J(q-\epsilon d)}{2\epsilon}
$$

与伴随梯度方向导数：

$$
g^T d
$$

比较。

代码为：

```matlab
finiteDifference = (plusObjective-minusObjective)/(2*epsilon);
adjointDirection = gradient'*direction;
```

### 19.3 Gauss--Newton 二次型检查

比较观测增量的平方：

```matlab
tangentObservation'*tangentObservation
```

与有限差分观测增量的平方：

```matlab
finiteDifferenceObservation'*finiteDifferenceObservation
```

这验证的是 Gauss--Newton Hessian 中的观测 Jacobian 部分。

## 20. 画图部分

脚本最后生成若干图。

### Figure 1: beta 恢复结果

画出：

- `betaTrue`；
- `betaInitial`；
- `betaRecovered`；
- 参数节点上的恢复值。

对应代码：

```matlab
plot(plotFine,betaTruePlot,'k-',...)
plot(plotFine,betaInitialPlot,'b--',...)
plot(plotFine,betaRecoveredPlot,'r-',...)
plot(xBeta,betaRecovered,'ro',...)
```

### Figure 2: beta 误差历史

画出：

- 相对二范数误差；
- 绝对无穷范数误差；
- 相对无穷范数误差。

### Figure 3: 顶部速度观测拟合

比较：

- synthetic observation；
- recovered prediction。

也就是：

```matlab
plot(xObs,dataObs,'ko',xObs,uRecovered(topDof),'r-',...)
```

### Figure 4: 恢复解场

画出：

- 恢复的水平速度 $u_x$；
- 恢复的竖直速度 $u_y$；
- 速度大小 $|u|$；
- 压力 $p$。

这些图用于检查恢复后的正问题解是否合理。

## 21. 局部函数说明

### `solveforward`

```matlab
function [u,eqn,info,p] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L)
```

作用：

1. 把 `q` 转换成 `beta=exp(q)`；
2. 构造周期 P1 插值函数；
3. 如果有 warm start，则放入 `option.u0`；
4. 调用 `NonlinearStokesP2P1`；
5. 返回速度、压力、方程结构和求解信息。

### `assembleparameterderivative`

构造

$$
G=R_q.
$$

它逐个参数方向调用：

```matlab
eqn.applyBetaDerivative(directionFunction)
```

得到残差对参数方向的导数。

### `gaussnewtonproduct`

给定一个参数方向 `direction`，返回

$$
\left(H_{\mathrm{GN}}+\lambda I\right)\,\mathrm{direction}.
$$

这是 PCG 所需的矩阵自由乘法函数。

### `verifyderivatives`

用有限差分检查：

- 增量状态方程；
- 伴随梯度；
- Gauss--Newton 观测 Jacobian。

### `periodicP1`

实现周期线性插值。用于把参数节点值变成空间函数。

### `boundaryweights`

构造顶部边界积分权重。

### `trimhistory`

如果反演提前停止，裁剪历史数组，只保留已经执行的迭代。

## 22. 和其他脚本的关系

### `NSDiagnosis.m`

这个脚本用于诊断不同 Glen 指数下的反演问题，并使用有限差分 Jacobian 做 SVD 分析。它更适合理解问题病态性和观测可辨识性。

相比之下，`NonlinearStokesAdjointInversionBoundaryObjective.m` 使用伴随梯度和矩阵自由 Gauss--Newton，更接近实际反演算法。

### `NSDerivativeComparison.m`

这个脚本用于验证：

- FD 梯度与伴随梯度是否一致；
- FD 构造的 Gauss--Newton Hessian 与 tangent/adjoint 构造是否一致。

它是对 `NonlinearStokesAdjointInversionBoundaryObjective.m` 中导数实现的独立检查。

### `NonlinearStokesP2P1.m`

这是正问题求解器。它提供：

- 非线性 Picard 求解；
- 一致切线矩阵 `eqn.tangent`；
- 参数导数接口 `eqn.applyBetaDerivative`。

反演脚本依赖它完成正问题和导数相关计算。

## 23. 常见问题

### `eqn.tangent` 是什么？

它是完整非线性离散残差对状态变量的导数：

$$
eqn.tangent = R_U(U,q).
$$

这里的状态变量包括速度、压力和约束乘子。

它不同于 Picard 迭代中冻结黏度得到的矩阵。`eqn.tangent` 是一致切线矩阵，包含黏度和滑移系数对速度的导数。

### 为什么用 `q=log(beta)`？

因为底部滑移系数必须为正。用

$$
\beta=\exp(q)
$$

后，无论 `q` 如何更新，都有 $\beta>0$。

### 为什么 Hessian 是 Gauss--Newton Hessian？

目标函数是观测残差平方：

$$
J(q)=\frac{1}{2}\|r(q)\|^2.
$$

完整 Hessian 为

$$
\nabla^2J
=
J_r^T J_r
+
\sum_i r_i\nabla^2 r_i.
$$

Gauss--Newton 方法只保留第一项：

$$
H_{\mathrm{GN}}=J_r^T J_r.
$$

脚本中的 `gaussnewtonproduct` 实现的是这一近似，而不是完整 Hessian。

### 为什么 PCG 可以用函数句柄？

PCG 只需要计算矩阵作用在向量上的结果，不需要显式矩阵。

因此代码写成：

```matlab
hessian = @(direction) gaussnewtonproduct(...);
pcg(hessian,-gradient,...)
```

这样可以避免显式装配 Hessian。

### `fPicard` 和 `pcgIt` 有什么区别？

`fPicard` 是非线性正问题 Picard 迭代步数。

`pcgIt` 是优化步中解线性 Gauss--Newton 系统的 PCG 迭代步数。

它们属于不同层次：

- Picard：求解非线性 PDE；
- PCG：求解优化子问题中的线性系统。

### 为什么 objective 很小后 beta 误差不一定为零？

顶部速度观测不一定能完全区分所有底部参数模式。某些 $\beta$ 的变化对顶部速度影响很弱，导致数据拟合已经很好，但参数仍有误差。

这也是 `NSDiagnosis.m` 做 SVD 诊断的原因。

## 24. 一句话总结

`NonlinearStokesAdjointInversionBoundaryObjective.m` 是一个基于顶部边界速度积分误差的非线性 Stokes 底部滑移反演脚本。它用 Picard 方法求正问题，用一致切线矩阵解伴随方程得到梯度，用矩阵自由 Gauss--Newton/LM 方法计算更新步，并通过 line search 保证目标函数下降。
