# `NonlinearStokesAdjInvSlabBed.m` 说明

本文解释脚本

```text
ice-sheet/FullStokes2d/NonlinearStokesAdjInvSlabBed.m
```

的数学目标、代码流程和几个关键导数对象。这个脚本做的是二维非线性 Stokes 冰流模型中的底部滑移系数反演：

- 正问题：给定底部滑移系数 $\beta(x)$，求速度 $u$ 和压力 $p$；
- 观测：只使用顶部边界上的水平速度；
- 反问题：从顶部速度观测恢复底部滑移系数；
- 优化变量：不是直接优化 $\beta$，而是优化

$$
q=\log(\beta),
\qquad
\beta=\exp(q).
$$

这样可以保证 $\beta>0$。脚本没有加入额外正则项，核心目标函数是顶部边界速度误差：

$$
J(q)
=
\frac12
\frac{
\int_{\Gamma_t}(u(q)-u_{\mathrm{obs}})^2\,\mathrm{d}s
}{
\int_{\Gamma_t}u_{\mathrm{obs}}^2\,\mathrm{d}s
}.
$$

整体算法可以概括为：

1. 用真实参数 `qTrue` 解一次正问题，生成合成观测；
2. 从扰动初值 `q` 出发反演；
3. 每轮用阻尼 Picard 迭代解非线性 Stokes 正问题；
4. 用一致切线矩阵 `eqn.tangent` 解伴随方程，得到梯度；
5. 用矩阵自由 Gauss--Newton / Levenberg--Marquardt 方法求更新步；
6. 更新参数并输出误差、收敛信息和图像。

## 1. 几何、边界和观测

脚本先在矩形区域上生成网格：

```matlab
L = 4;
H = 1;
h = 0.1;
slope = 0.1;

[node,elem] = squaremesh([0,L,0,H],h);
```

边界条件设置为：

```matlab
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
```

含义是顶部为 Neumann 边界，底部为 Robin 滑移边界。随后做坐标变换：

```matlab
node(:,2) = node(:,2)-slope*node(:,1);
```

所以实际区域是一个倾斜 slab，顶部边界为

$$
y=H-\mathrm{slope}\,x.
$$

脚本还启用周期边界：

```matlab
option.periodic = true;
option.periodic_x = [0,L];
```

速度使用 P2 元，因此速度自由度包含网格顶点和边中点：

```matlab
[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
```

顶部观测自由度 `topDof` 是位于顶部边界上的 P2 速度自由度。右端点被排除，因为周期条件会把 $x=L$ 和 $x=0$ 识别为同一个位置：

```matlab
topDof = find(abs(uNode(:,2)-surfaceLevel)<tolGeometry ...
            & uNode(:,1)<L-tolGeometry);
```

顶部边界积分用简单等权重近似：

```matlab
topWeight = boundaryweights(xObs,L,slope);
```

由于顶部是倾斜直线，弧长微元为

$$
\mathrm{d}s=\sqrt{1+\mathrm{slope}^2}\,\mathrm{d}x,
$$

所以每个观测点的权重近似为

$$
w_i=\sqrt{1+\mathrm{slope}^2}\frac{L}{N_{\mathrm{obs}}}.
$$

代码中

```matlab
topWeight'*(residual.^2)
```

就是对顶部边界积分的离散近似。

## 2. 正问题和参数化

非线性 Stokes 模型参数为：

```matlab
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];
```

其中 `pde.n` 是 Glen 指数，`pde.m` 是 Weertman 滑移律指数。正问题由

```matlab
NonlinearStokesP2P1(node,elem,bdFlag,pde,option)
```

求解，使用 P2 速度、P1 压力和阻尼 Picard 迭代。反演需要导数，所以设置：

```matlab
option.assemble_tangent = true;
```

底部滑移系数使用周期 P1 参数化：

```matlab
Nm = round(L/h);
xBeta = (0:Nm-1)'*h;
beta = exp(q(:));
pde.beta = @(pt) periodicP1(pt(:,1),xBeta,beta,L);
```

当前默认 `L = 4`、`h = 0.1`，所以

```text
Nm = 40
```

也就是有 40 个底部滑移参数自由度。`periodicP1` 使用线性插值，并令 $x=L$ 处的值等于 $x=0$ 处的值。

脚本先用真实参数 `qTrue` 生成合成观测：

```matlab
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
dataObs = uTrue(topDof);
```

这是一种 inverse-crime 风格的测试：观测数据和反演使用同一张网格，适合检查导数、伴随和优化流程是否一致。

## 3. 目标函数和每轮迭代

每轮反演迭代先解当前参数下的正问题：

```matlab
[u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
    node,elem,bdFlag,xBeta,L);
```

`uWarm` 是上一轮接受的速度解，用作 Picard 迭代初值。顶部观测残差为：

```matlab
residual = u(topDof)-dataObs;
```

离散目标函数为：

```matlab
dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
objective = dataObjective;
```

也就是

$$
J(q)
=
\frac12
\frac{\sum_i w_i(u_i-u_{\mathrm{obs},i})^2}
{\sum_i w_i u_{\mathrm{obs},i}^2}.
$$

其中

```matlab
dataNormSquared = max(topWeight'*(dataObs.^2),eps);
```

用于归一化，避免目标函数受速度量级影响。

接下来脚本计算参数导数矩阵 `G`，构造目标函数对状态变量的导数，解伴随方程，得到梯度：

```matlab
G = assembleparameterderivative(eqn,q,xBeta,L,Nm);

observationGradient = zeros(size(eqn.tangent,1),1);
observationGradient(topDof) = topWeight.*residual/dataNormSquared;

adjoint = eqn.tangent'\(-observationGradient);
gradient = G'*adjoint;
```

这几行是脚本的核心。下面分别解释 `R`、`eqn.tangent`、`G` 和 `applyBetaDerivative`。

## 4. 残差 `R` 和一致切线矩阵 `eqn.tangent`

把离散正问题写成

$$
R(U,q)=0.
$$

这里 `R` 是完整离散方程残差，也就是把当前状态代入方程后，各个离散方程还差多少才能平衡。它不是有限元解误差，而是离散方程的不平衡量。状态变量为

$$
U=
\begin{bmatrix}
u_x\\
u_y\\
p\\
\lambda
\end{bmatrix},
$$

包含二维速度、压力和约束乘子。残差可以理解为

$$
R(U,q)
=
\begin{bmatrix}
R_{\mathrm{momentum}}\\
R_{\mathrm{div}}\\
R_{\mathrm{constraint}}
\end{bmatrix},
$$

分别对应动量方程、不可压缩条件和周期/法向/压力约束。

`eqn.tangent` 是残差对状态变量的 Jacobian：

$$
eqn.tangent
=
R_U(U,q)
=
\frac{\partial R}{\partial U}.
$$

它描述状态变量发生微小变化时，完整残差的一阶变化：

$$
R(U+\delta U,q)
\approx
R(U,q)+R_U\delta U.
$$

令 `Nu` 为一个速度分量的 P2 自由度个数，`Np` 为 P1 压力自由度个数，`nConstraint` 为约束个数，则

$$
eqn.tangent
\in
\mathbb{R}^{(2Nu+Np+nConstraint)\times(2Nu+Np+nConstraint)}.
$$

它由 `NonlinearStokesP2P1.m` 在正问题收敛后装配。相关代码逻辑是：

```matlab
if option.assemble_tangent
    Kt = assembleviscoustangent(u);
    Kbt = assemblebedtangent(u);
    tangentM = [Kt+Kbt, B'; B, sparse(Np,Np)];
    eqn.tangent = [tangentM, C'; C, sparse(nConstraint,nConstraint)];
    eqn.applyBetaDerivative = @assemblebetadirection;
end
```

对应分块结构为

$$
R_U
=
\begin{bmatrix}
K_{\mathrm{visc}}^{\mathrm{tan}}
+K_{\mathrm{bed}}^{\mathrm{tan}} & B^T & C^T\\
B & 0 & 0\\
C & 0 & 0
\end{bmatrix}.
$$

其中：

- `Kt = assembleviscoustangent(u)` 是 Glen 黏性残差的一致切线；
- `Kbt = assemblebedtangent(u)` 是底部滑移残差的一致切线；
- `B` 是速度-压力耦合块；
- `C` 是周期条件、底部法向约束、压力规范化等约束块。

`eqn.tangent` 不是 Picard 迭代中的冻结系数矩阵。Picard 矩阵适合稳健求解正问题，但通常不是完整非线性残差的 Jacobian。非线性黏性项为

$$
r_{\mathrm{visc}}(u;v)
=
\int_\Omega
2\eta(u)\dot{\varepsilon}(u):\dot{\varepsilon}(v)\,\mathrm{d}x.
$$

对速度方向 $\delta u$ 求导时，有两部分：

$$
\begin{aligned}
D r_{\mathrm{visc}}(u)[\delta u;v]
= {}&
\int_\Omega
2\eta(u)
\dot{\varepsilon}(\delta u):\dot{\varepsilon}(v)\,\mathrm{d}x\\
&+
\int_\Omega
2D\eta(u)[\delta u]
\dot{\varepsilon}(u):\dot{\varepsilon}(v)\,\mathrm{d}x.
\end{aligned}
$$

Picard 矩阵通常冻结 $\eta(u)$，会漏掉第二项；一致切线矩阵必须包含它。底部非线性滑移也类似，如果牵引随切向速度变化，`Kbt` 必须包含滑移牵引对速度的导数。伴随梯度和灵敏度方程需要的是这个真实 Jacobian，所以脚本使用 `eqn.tangent`。

## 5. 参数导数矩阵 `G` 和 `applyBetaDerivative`

`G` 是残差对反演参数 `q` 的 Jacobian：

$$
G
=
R_q(U,q)
=
\frac{\partial R}{\partial q}.
$$

它的维数是

$$
G
\in
\mathbb{R}^{(2Nu+Np+nConstraint)\times Nm}.
$$

当前默认 `Nm = 40`，所以 `G` 有 40 列。第 `j` 列表示只扰动第 `j` 个参数自由度 `q_j` 时，完整残差的一阶变化。

`assembleparameterderivative` 逐列装配 `G`：

```matlab
function G = assembleparameterderivative(eqn,q,xBeta,L,Nm)
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        deltaBeta = beta.*direction;
        directionFunction = @(pt) periodicP1(...
            pt(:,1),xBeta,deltaBeta,L);
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end
```

关键是链式法则。正问题依赖的是 $\beta$，优化变量是 $q$，且

$$
\beta=\exp(q).
$$

因此

$$
\delta\beta
=
\operatorname{diag}(\beta)\delta q.
$$

当 $\delta q=e_j$ 时，只有第 `j` 个参数节点有扰动，大小为 $\beta_j$。代码中的

```matlab
deltaBeta = beta.*direction;
```

正是这个转换。

`directionFunction` 把离散的 `deltaBeta` 节点值转换成可在底部积分点求值的周期 P1 函数。随后

```matlab
G(:,j) = eqn.applyBetaDerivative(directionFunction);
```

计算这一方向的残差导数。

`eqn.applyBetaDerivative` 是 `NonlinearStokesP2P1.m` 返回的函数句柄，指向局部装配函数 `assemblebetadirection`。它做的是：

$$
\delta\beta(x)
\longmapsto
R_\beta\,\delta\beta.
$$

也就是说，它不是重新求解正问题，而是在当前收敛状态 $U$ 上，计算一个给定底部滑移扰动 $\delta\beta$ 对残差的一阶影响。

底部牵引可写成

$$
t_b(u_t,\beta)
=
\beta s_b^{(m-1)/2}u_t.
$$

保持当前速度 $u$ 不变，只扰动 $\beta$，有

$$
D_\beta t_b[\delta\beta]
=
\delta\beta\,s_b^{(m-1)/2}u_t.
$$

把这个量乘以有限元测试函数并在底部边界积分，就得到 `applyBetaDerivative` 返回的向量。二者关系可以概括为：

```text
applyBetaDerivative:         delta beta(x) -> R_beta delta beta
assembleparameterderivative: 所有 q_j 方向 -> G = R_q
```

因为优化变量是 `q`，所以最终使用的是

$$
R_q\delta q
=
R_\beta(\beta\,\delta q).
$$

## 6. 伴随梯度

目标函数只依赖顶部观测速度，所以 `J_U` 只有 `topDof` 上非零：

```matlab
observationGradient = zeros(size(eqn.tangent,1),1);
observationGradient(topDof) = topWeight.*residual/dataNormSquared;
```

伴随方程为

$$
R_U^T z=-J_U^T.
$$

代码是：

```matlab
adjoint = eqn.tangent'\(-observationGradient);
```

然后 reduced gradient 为

$$
\nabla J(q)=G^T z.
$$

代码是：

```matlab
gradient = G'*adjoint;
```

维数上，

$$
z\in\mathbb{R}^{2Nu+Np+nConstraint},
\qquad
G^Tz\in\mathbb{R}^{Nm}.
$$

所以 `gradient` 和参数向量 `q` 维数一致。

脚本没有再做 Picard 迭代来求伴随，因为伴随问题是在收敛后的状态上解一个线性系统。它不是重新求一个非线性正问题，而是直接使用已有的 `eqn.tangent'`。

## 7. Gauss--Newton / LM 更新

脚本用 PCG 解阻尼 Gauss--Newton 方程：

$$
(H_{\mathrm{GN}}+\lambda I)s=-g.
$$

这里 `lambda` 是 Levenberg--Marquardt 阻尼参数。PCG 不需要显式矩阵，只需要 Hessian-vector product：

```matlab
hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
    topDof,topWeight,dataNormSquared,lambda);

[step,flag,relativeResidual,pcgIt] = pcg(...
    hessian,-gradient,pcgTolerance,pcgMaxIt);
```

`gaussnewtonproduct` 给定一个参数方向

$$
direction=\delta q\in\mathbb{R}^{Nm},
$$

先解增量状态方程：

```matlab
incrementalState = eqn.tangent\(-G*direction);
```

即

$$
R_U\delta U=-R_q\delta q.
$$

再把状态增量投影到观测空间：

```matlab
incrementalObservation(topDof) = ...
    topWeight.*incrementalState(topDof)/dataNormSquared;
```

然后解增量伴随：

```matlab
incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
```

最后返回

```matlab
product = G'*incrementalAdjoint+lambda*direction;
```

也就是参数空间中的

$$
(H_{\mathrm{GN}}+\lambda I)\delta q.
$$

这套写法的好处是不用显式形成稠密的观测 Jacobian，也不用显式形成 Hessian。

## 8. 步长、停止条件和输出

默认情况下

```matlab
useLineSearch = defaultlinesearch();
```

而 `defaultlinesearch()` 返回 `false`，所以脚本直接尝试完整 Gauss--Newton / LM 步：

```matlab
qTrial = q+step;
```

如果 trial 正问题收敛，就接受该步，并减小阻尼：

```matlab
lambda = max(lambda/3,1e-12);
```

如果 trial 正问题不收敛，则不更新 `q`，并增大阻尼：

```matlab
lambda = 10*lambda;
```

若手动把 `defaultlinesearch()` 改成 `true`，脚本会做最多 10 次回溯线搜索，要求 trial objective 下降才接受。

主要停止条件包括：

- `objective < 1e-15`：目标函数已经非常小；
- `norm(gradient) <= gradientTolerance`：梯度足够小；
- `norm(step) <= stepTolerance*max(1,norm(q))`：步长足够小；
- 达到 `maxInverseIt`。

迭代表格字段含义：

- `it`：反演迭代编号；
- `objective`：当前目标函数；
- `betaLinfRel`：$\beta$ 的相对 $L^\infty$ 误差；
- `|grad|`：参数梯度二范数；
- `fPicard`：本轮正问题 Picard 迭代步数；
- `pcgIt`：PCG 迭代步数；
- `pcgRel`：PCG 返回的相对残差；
- `ls`：线搜索次数；
- `time`：本轮耗时；
- `stop`：停止原因。

## 9. 导数检查和画图

第一轮迭代会调用

```matlab
verifyderivatives(...)
```

做三个有限差分检查：

1. 增量状态方程是否正确：

   $$
   \delta U
   =
   -R_U^{-1}R_q\delta q;
   $$

2. 伴随梯度方向导数是否和有限差分一致：

   $$
   \nabla J(q)^T\delta q
   \approx
   \frac{J(q+\varepsilon\delta q)-J(q-\varepsilon\delta q)}
   {2\varepsilon};
   $$

3. Gauss--Newton 观测二次型是否和有限差分观测增量一致。

这些检查用于确认 `eqn.tangent`、`G`、伴随梯度和 Gauss--Newton product 是相互一致的。

脚本最后绘制四类图：

- Figure 1：真实、初始和恢复的 $\beta$；
- Figure 2：目标函数和 $\beta$ 误差历史；
- Figure 3：顶部速度观测与恢复预测；
- Figure 4：网格、恢复的 $u_x$、$u_y$ 和压力。

## 10. 局部函数概览

`solveforward` 把 `q` 转换成 `beta = exp(q)`，构造周期 P1 插值函数，设置 warm start，然后调用 `NonlinearStokesP2P1`。

`assembleparameterderivative` 逐列构造

$$
G=R_q.
$$

它通过链式法则 $\delta\beta=\beta\,\delta q$，把每个 `q_j` 方向转换成对应的底部滑移扰动，再调用 `eqn.applyBetaDerivative`。

`gaussnewtonproduct` 实现矩阵自由的

$$
(H_{\mathrm{GN}}+\lambda I)\delta q.
$$

`verifyderivatives` 用有限差分检查增量状态、伴随梯度和 Gauss--Newton 二次型。

`periodicP1` 负责把周期参数节点值插值到任意积分点。

`boundaryweights` 构造顶部边界积分权重。

`trimhistory` 在提前停止时裁剪历史数组。

## 11. 一句话总结

`NonlinearStokesAdjInvSlabBed.m` 用顶部速度误差反演底部滑移参数。正问题用 Picard 迭代求解，但反演导数必须使用完整非线性残差的 Jacobian：

```text
eqn.tangent = R_U
G           = R_q
```

伴随梯度由

```matlab
adjoint = eqn.tangent'\(-observationGradient);
gradient = G'*adjoint;
```

给出，Gauss--Newton 更新则通过增量状态和增量伴随方程以矩阵自由方式计算。
