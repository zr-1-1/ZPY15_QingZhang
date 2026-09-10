# Data_Extract 轨道递推耗时与优化方案

日期：2026-09-10。对象：[Data_Extract.m](../src/Data_Extract.m) 第 177～192 行“递推所有碎片轨道，并保存，方便后续导入”，以及 [OE_scl_ptb.m](../src/OE_scl_ptb.m) 和脚本内的轨道转换函数。

本文记录本机抽样基准结果与优化建议。优化算法在临时脚本中验证，尚未应用到项目源码。计时不包括 ATK 数据提取、可视化和目标节之后的计算。

## 1. 结论

**当前代码预计需要 2.5～3 分钟；优先采用“逐碎片处理、时间向量化”，预计可将计算与保存合计降至 20～30 秒。多线程适合作为后续可选优化。**

| 方案 | 20 个碎片、完整 24 小时实测 | 345 个碎片计算耗时外推 |
|---|---:|---:|
| 当前逐点双层 `for` | 8.065 秒 | 139.13 秒 |
| 时间向量化，并提取固定参数、展开旋转矩阵 | 0.117～0.137 秒 | 2.01～2.37 秒 |
| 上述向量化 + 4 线程，预热后 | 0.080 秒 | 1.38 秒 |

当前保存方式的样本耗时为 0.991 秒，按碎片数量线性外推，完整保存约 17.10 秒。向量化后，保存将占据大部分时间。

这些是抽样外推，未执行 345 个碎片的完整优化任务。全量数组的分配、内存访问、磁盘负载和 MATLAB 当前负载都会影响实际结果；表内小数仅保留测量依据，不代表预测具有同等精度。

## 2. 当前工作量与存储规模

当前参数为 345 个碎片、24 小时、1 秒采样一次：

```matlab
t_step = 1;
Tr_0 = 60*60*24;
x_rv = zeros(6,345,Tr_0/t_step+1);
```

- 待计算状态数：`345 × 86400 = 29,808,000`。
- 包含初始帧的输出尺寸：`6 × 345 × 86401`。
- 输出为 `double`，数组净数据量：`6 × 345 × 86401 × 8 = 1,430,800,560` 字节，约 **1.431 GB / 1.333 GiB**。
- 每个时间点调用一次 `OE_scl_ptb` 和一次 `orbit6toRV`，后者包含开普勒方程迭代和坐标转换。

数组净数据量不等于 MATLAB 峰值内存。中间数组、并行结果与保存缓冲都可能增加占用。

当前传播使用初始根数 `E_0` 与累计时间 `t` 直接计算各时刻状态，各时刻相互独立。它不是逐秒调用数值积分器，也不依赖上一个时间点的结果。

## 3. 基准方法与结果

### 3.1 环境与输入

- MATLAB：`24.2.0.2712019 (R2024b)`，`PCWIN64`。
- 输入：项目 `data/debris_orbits.csv`，共 345 行。
- 根数初始化：从当前脚本提取 `rv2orb_elements`，使用 CSV 的位置与速度计算。
- 参数：`mu = 3.986004418e14`、`R_E = 6378137`、`J_2 = -sqrt(5)*(-4.841653717360e-04)`。
- 抽样索引：`round(linspace(1,345,20))`。这是按行号均匀抽样，不是按轨道参数分层抽样。
- 时间点：`1:86400`，每种方案计算 1,728,000 个状态。
- 在独立 MATLAB 批处理进程中运行，保存样本写入临时目录。

基线使用当前标量传播与转换函数，保留逐点赋值方式，样本数组为 `6×20×86401`。数组分配和输入准备位于计算计时之外，省略原代码每个碎片一次的 `i` 显示。样本第 0 秒帧保留为零，误差比较仅覆盖第 1～86400 秒；实际实现必须保留 CSV 初始状态。

### 3.2 测量记录

| 项目 | 实测结果 |
|---|---:|
| 标量基线，较早一次独立运行 | 7.952633 秒 |
| 标量基线，最终对照运行 | 8.065225 秒 |
| 默认 `save` 保存最终基线样本 | 0.991179 秒 |
| 样本 MAT 文件大小 | 约 80.400 MB |
| 向量化，第 1 次 | 0.137273 秒 |
| 向量化，第 2 次 | 0.118285 秒 |
| 向量化，第 3 次 | 0.116748 秒 |
| 4 线程池启动 | 1.325 秒 |
| 4 线程向量化，第 1 次执行 | 2.007434 秒 |
| 4 线程向量化，第 2 次执行 | 0.080010 秒 |

外推方法：`样本耗时 × 345/20`。保存文件大小同样外推约为 **1.387 GB**，实际大小取决于数据和保存设置。

并行首次执行含有初始化、调度等固定开销，不能将首次耗时简单乘以 `345/20`。表中并行计算外推仅采用第二次执行，且不包括线程池启动。当前串行向量化已很快，单次任务未必能通过并行缩短总时间。

### 3.3 数值一致性

向量化版本保持原开普勒求解阈值 `0.001*pi/180`，通过活动掩码让每个时间点独立停止牛顿迭代，保留原来的真近点角转换和 `E == pi` 分支。

对 20 个样本的全部 86,400 个计算时刻，与标量基线逐分量比较：

| 检查项 | 最大绝对差 |
|---|---:|
| 位置分量 | 0 m |
| 速度分量 | `1.81898940355e-12 m/s` |

这是优化前后的数值一致性验证，不是对 J₂ 平均摄动模型的物理精度验证。误差检查针对串行向量化结果；线程版本复用了相同计算函数，本次未另做线程输出误差检查。

## 4. 优化方案与顺序

### 4.1 首选：按时间向量化

保留外层碎片循环，每次批量计算一个碎片的所有时间点：

1. 生成时间向量 `time_s = 0:t_step:Tr_0`。
2. 从该碎片的根数中提取常量，批量求出升交点赤经、近地点幅角和平近点角。
3. 对时间向量批量求解开普勒方程。
4. 批量转换为位置速度，并一次写入该碎片的输出切片。

开普勒求解不能让已收敛时间点持续迭代，否则可能与原标量停止规则产生差异。应保留逐点收敛状态；实际实现还应对非法输入和不收敛情形给出明确处理。

### 4.2 同时提取常量并简化坐标转换

当前代码重复计算平均运动 `n`、摄动系数 `C_J2`、偏心率相关分母和倾角三角函数。这些在同一个碎片的整段传播中不变，应移至时间计算之外。

`orbit6toRV` 每次构造三个 `3×3` 旋转矩阵再相乘，可将位置变换展开为逐元素公式，避免数千万次小矩阵构造和乘法。

本次向量化提速包含上述改动，不能将全部提速单独归因于移除内层 `for`。

建议先逐碎片向量化，避免同时生成所有碎片、所有时间点的多份大型中间数组。单个碎片的完整六维 `double` 轨迹约为 4.15 MB，更容易控制临时内存。

### 4.3 可选：外层 `parfor` 与线程池

按碎片分配并行任务；每个任务内部仍然采用时间向量化。示意接口如下，`propagate_vectorized` 是待实现函数，不是现有项目接口：

```matlab
% oe: 6×N；Debris_rv: 6×N。
% time_s 须经过步长、末时刻等约束校验。
time_s = 0:t_step:Tr_0;
N = size(oe,2);
x_rv = zeros(6,N,numel(time_s));

% 调用前按需建立或复用线程池，例如 parpool('Threads',4)。
parfor i = 1:N
    trajectory = zeros(6,numel(time_s));
    trajectory(:,1) = Debris_rv(:,i);
    trajectory(:,2:end) = propagate_vectorized( ...
        oe(:,i),time_s(2:end),J_2,R_E,mu);
    x_rv(:,i,:) = reshape(trajectory,6,1,[]);
end
```

实施要点：

- 将 `J_2`、`R_E`、`mu` 显式传入工作函数，避免依赖工作线程中的 `global` 初始化状态。
- 使用完整的 `x_rv(:,i,:)` 切片写入，便于 `parfor` 分析输出变量。
- 将池启动与实际计算分别计时。若已有合适的池，优先复用，不在每次调用中重复创建和销毁。
- 本机已实际成功创建 4 线程池；没有比较更多线程数，也没有测试进程池。
- 没有测试“仅将原始外层 `for` 改为 `parfor`”的性能，不对该方案给出倍数保证。

**推荐默认采用串行向量化；反复批量运行时，再根据端到端实测决定是否启用线程池。**

### 4.4 后续重点：保存耗时

向量化后，约 17 秒的保存外推时间远高于约 2 秒的计算时间，应把进一步优化重点放到输出阶段。

可以用同一完整数据分别测试当前 `save` 与以下形式，记录耗时、文件大小及读回一致性：

```matlab
save(outputFile,'x_rv','-nocompression');
```

关闭压缩的收益尚未实测，不纳入本文提速结论。存储速度和压缩开销会共同决定结果。若后续有局部读取或更大变量需求，再单独评估 `-v7.3`、`matfile` 和分块方案，不能预设它们比当前整块保存更快。

## 5. 实施时的验收要求

性能优化应保留当前输出接口和物理含义：

- 输出仍为 `6×N×Nt`，前 3 行位置、后 3 行速度，单位分别为 m 和 m/s。
- 初始帧直接使用 `Debris_rv`；碎片顺序与输入一致。
- 时间与存储索引分离，支持非 1 秒步长时，不能继续使用 `t+1` 作为下标。
- 明确根数矩阵布局，避免依赖重复执行节首转置来切换方向。
- 保持原模型、单位和开普勒收敛规则；任何额外精度调整单独验证。
- 完整运行后检查全部结果有限，分别比较位置和速度误差，并验证保存读回一致性。
- 分别记录初始化、计算、并行池启动和保存时间；最终以完整 345 个碎片的端到端耗时验收。

提高 `t_step` 虽然减少计算量和文件大小，但会改变采样分辨率，不作为保持当前输出要求的首选优化。改用 `single` 也会改变精度，本次未采用。

## 6. 可直接使用的优化代码

下面代码分为目标节替换代码和独立函数文件。默认保持当前压缩保存方式，采用串行向量化；将 `useThreads` 改为 `true` 即可使用 4 线程。代码仍使用 `double` 和原 J₂ 平均摄动模型。

### 6.1 替换原“递推所有碎片轨道，并保存，方便后续导入”节

需要先完成路径、全局常量和 CSV 初始状态的初始化。将下一小节的函数保存为 `src/propagate_debris_history.m`，确保 `src` 已加入 MATLAB 路径，再使用本段：

```matlab
%% 递推所有碎片轨道，并保存，方便后续导入
t_step = 1;
Tr_0 = 24*60*60;
useThreads = false;  % 首选串行向量化；反复批量运行时可设为 true

% 函数接受 N×6 或 6×N 根数，不再原地转置 Debris_oe。
% 池启动包含在本段计时中，另有独立的启动耗时输出。
computeTimer = tic;
[x_rv,time_s] = propagate_debris_history( ...
    Debris_oe,Debris_rv,Tr_0,t_step,J_2,R_E,mu,useThreads);
fprintf('轨迹生成（含初始化及按需启动并行池）：%.3f s\n',toc(computeTimer));

saveTimer = tic;
save(fullfile(projectRoot,'data','x_rv.mat'),'x_rv');
fprintf('保存：%.3f s\n',toc(saveTimer));
% 如后续导入也需要时间轴，可将 time_s 加入 save 的变量列表。
```

本节不改变 `Debris_oe` 本身的方向。相邻后续代码如果依赖旧节的转置副作用，需要显式统一其输入布局。

### 6.2 保存为 `src/propagate_debris_history.m`

```matlab
function [x_rv,time_s] = propagate_debris_history( ...
    oe,rv0,duration_s,step_s,J2,RE,mu,useThreads)
% 输出 x_rv: 6×N×Nt；位置 m，速度 m/s；time_s 从 0 开始。
% oe: N×6 或 6×N，根数顺序 [a,e,i,Omega,omega,M]，单位 m/rad。
% rv0: 6×N，直接作为初始帧。6×6 的 oe 按每列一个碎片解释。
if nargin < 8
    useThreads = false;
end
validateattributes(rv0,{'double'},{'real','finite','nonempty','2d'});
assert(size(rv0,1)==6,'rv0 必须为 6×N。');
N = size(rv0,2);
validateattributes(oe,{'double'},{'real','finite','nonempty','2d'});
if isequal(size(oe),[6,N])
    % 已是每列一个碎片。
elseif isequal(size(oe),[N,6])
    oe = oe.';
else
    error('oe 尺寸必须为 6×N 或 N×6，并与 rv0 的碎片数一致。');
end
assert(all(oe(1,:)>0) && all(oe(2,:)>=0 & oe(2,:)<1), ...
    '仅支持 a>0、0<=e<1 的椭圆轨道。');
validateattributes(duration_s,{'double'},{'real','finite','scalar','nonnegative'});
validateattributes(step_s,{'double'},{'real','finite','scalar','positive'});
validateattributes(J2,{'double'},{'real','finite','scalar'});
validateattributes(RE,{'double'},{'real','finite','scalar','positive'});
validateattributes(mu,{'double'},{'real','finite','scalar','positive'});
validateattributes(useThreads,{'logical'},{'scalar'});
steps = duration_s/step_s;
assert(isfinite(steps) && steps<=flintmax && ...
    abs(steps-round(steps))<=8*eps(max(1,steps)), ...
    '总时长必须是步长的整数倍。');
Nt = round(steps)+1;
time_s = (0:Nt-1)*step_s;
x_rv = zeros(6,N,Nt);
if Nt==1
    x_rv(:,:,1) = rv0;
    return
end
t = time_s(2:end);

if useThreads
    pool = gcp('nocreate');
    if isempty(pool)
        poolTimer = tic;
        pool = parpool('Threads',4);
        fprintf('线程池启动：%.3f s\n',toc(poolTimer));
    end
    assert(isa(pool,'parallel.ThreadPool'), ...
        '已有非线程并行池；请使用串行模式，或自行切换为线程池。');
    parfor k = 1:N
        trajectory = [rv0(:,k),propagate_one(oe(:,k),t,J2,RE,mu)];
        x_rv(:,k,:) = reshape(trajectory,6,1,Nt);
    end
else
    for k = 1:N
        trajectory = [rv0(:,k),propagate_one(oe(:,k),t,J2,RE,mu)];
        x_rv(:,k,:) = reshape(trajectory,6,1,Nt);
    end
end
end

function rv = propagate_one(q,t,J2,RE,mu)
% 每个碎片独立计算；固定系数仅求一次，时间维度向量化。
a = q(1); e = q(2); inc = q(3);
n = sqrt(mu/a^3);
c = 1.5*J2*(RE/a)^2*n;
O = q(4)-c*cos(inc)/(1-e^2)^2*t;
w = q(5)+c*(2-2.5*sin(inc)^2)/(1-e^2)^2*t;
M = q(6)+(n+0.5*c*(1-e^2)^(-1.5)*(3*cos(inc)^2-1))*t;
assert(all(isfinite(M)) && all(isfinite(O)) && all(isfinite(w)), ...
    '传播根数出现非有限值，请检查输入尺度。');

% 各时间点分别停止迭代，保持原 solve_Kepler 的停止规则。
% 新增 50 次上限，使异常输入明确报错，避免无限循环。
E = M;
active = true(size(E));
tol = 0.001*pi/180;
for iteration = 1:50
    old = E(active);
    next = old-(old-e*sin(old)-M(active))./(1-e*cos(old));
    assert(all(isfinite(next)),'开普勒迭代出现非有限值。');
    E(active) = next;
    active(active) = abs(next-old)>=tol;
    if ~any(active)
        break
    end
end
assert(~any(active),'开普勒方程在 50 次迭代内未收敛。');

theta = mod(2*atan(sqrt((1+e)/(1-e))*tan(E/2)),2*pi);
theta(E==pi) = pi;
distance = a*(1-e^2)./(1+e*cos(theta));
x = distance.*cos(theta);
y = distance.*sin(theta);
cw = cos(w); sw = sin(w);
co = cos(O); so = sin(O);
ci = cos(inc); si = sin(inc);

% 展开原 C3*C2*C1*r_o，避免逐时间点构造 3×3 矩阵。
rx = (co.*cw-so.*ci.*sw).*x+(-co.*sw-so.*ci.*cw).*y;
ry = (so.*cw+co.*ci.*sw).*x+(-so.*sw+co.*ci.*cw).*y;
rz = si*sw.*x+si*cw.*y;
u = w+theta;
A = sin(u)+e*sin(w);
B = cos(u)+e*cos(w);
speedScale = sqrt(mu/(a*(1-e^2)));
rv = [rx;ry;rz; ...
    speedScale*(-co.*A-so.*B*ci); ...
    speedScale*(-so.*A+co.*B*ci); ...
    speedScale*B*si];
assert(all(isfinite(rv),'all'),'输出轨迹出现非有限值。');
end
```

此封装增加了输入检查、初始帧拼接和迭代上限，故实际耗时应重新测量；第 3 节的提速数据对应临时算法基准，不是这份完整封装的端到端耗时。

## 7. 版本依据与适用范围

本次基准提取了当前源文件中的局部函数，并使用项目 CSV；未运行完整主脚本。先前的[检查报告](Data_Extract轨道递推与保存检查报告.md)和[复查报告](Data_Extract轨道递推与保存复查报告.md)对应各自检查时的版本，其中行号、布局问题及运行环境不能直接套用到当前版本。本文仅给出当前目标计算节的性能证据，不重新认定整个脚本的正确性。

本次记录的 SHA-256：

```text
src/Data_Extract.m
FCB4011FCC0F5906FBDA099D66A5ACE93F50B1ECABBEA2D49BC77BE437C0A4E1

src/OE_scl_ptb.m
C9EB2CABA47DE5FA9821C12EDFF950B1BB6C361A82853B02D5F7300DDF6C4C93

data/debris_orbits.csv
E2D9002F5C900D1568EC5EF5ADFC64427199A00680120F0CFA9A4DE26FF739D4
```
