# ZPY15_QingZhang
第十五届全国大学生周培源力学竞赛团体赛试题方案---清障小东风队针对清华方案的复现

> 2026-09-15：原 `data/mission` 的 113 个方案在 ATK 中出现高度违规，保留为历史结果。修正场景为 [`ZPY15_corrected.xml`](ZPY15_corrected.xml)，对应说明见 [ATK 修正与复核说明](docs/ATK修正与复核说明.md)。修正版离线校核为 112 个，尚待原生 ATK 重跑确认；不要继续将原表作为已验证的 ATK 方案。

# 当前仓库结构
```text
ZPY15_QingZhang/
├── ZPY15.xml                         # ATK 仿真场景文件
├── EGM96.grv                         # 重力场系数文件 EGM96，见赛题说明，位于"ATK-ZPY15专用版\AstroData\Earth\EGM96.grv"
├── docs/                             # 文档目录
└── src/
    ├── export_debris_csv.py          # 将 XML 中碎片初始状态导出为 CSV（Python 标准库）
    ├── Data_Extract.m                # 碎片轨道数据提取与计算脚本
    ├── OE_scl_ptb.m                  # J2 摄动下的轨道根数计算
    ├── orb_elements2rv.m             # 轨道根数转位置、速度矢量
    ├── rv2coe.m                      # 位置、速度矢量转轨道根数
    └── atk_connect_matlab_dependence/ # MATLAB 与 ATK 连接所需依赖
```
`src`存放脚本，`docs`存放文档

# 注意事项
1. 修改代码请先拉取远程仓库最新代码，避免冲突
2. .gitignore 文件用于忽略不需要上传的文件，如有需要请自行加入你不想上传的文件
3. 计算统一采用 SI 和弧度制：位置、半长轴和 `R_E` 为 m，速度为 m/s，时间为 s，`mu` 为 m³/s²，所有轨道角度为 rad，角速度为 rad/s，偏心率和 `J_2` 无量纲。

轨道根数统一按 `[a, e, i, Omega, omega, M]` 排列。`rv2coe` 和 `orb_elements2rv` 使用全局 `mu`；`OE_scl_ptb` 还使用全局 `J_2` 和 `R_E`，调用前应执行 `Data_Extract.m` 的参数初始化节。`R_E` 现在仅表示米制半径，不再使用 `R_E_m`。XML/CSV 的位置和速度保持 m、m/s；ATK 接口输入也须符合这一单位约定。已有 MAT 文件不自动换算，来源单位不明确时应重新提取。

`OE_scl_ptb` 返回累计弧度，允许角度超过 `2*pi`。可在 MATLAB 仓库根目录执行 `addpath('tests'); test_si_units` 运行单位与弧度验证。

# 常用命令
```bash
git clone <仓库地址>	    # 克隆一个远程仓库到本地，当前项目使用： git clone https://github.com/zr-1-1/ZPY15_QingZhang.git
git status	                # 查看工作区状态
# 同步上传相关命令
git add <文件名>	        # 将指定文件添加到暂存区
git add .	                # 添加当前目录下所有变更（新增、修改）到暂存区
git commit -m "提交说明"	# 将暂存区内容提交到本地仓库
git push origin main        # 将本地main分支推送到远程仓库 main 分支
git push <远程名> <分支名>	 # 将本地分支推送到远程仓库，远程名一般为origin
git push -u origin <分支名>	# 首次推送并建立本地分支与远程分支的关联
git pull <远程名> <分支名>	 # 拉取远程更新并自动合并到当前分支，远程名一般为origin，git pull 等同于 git fetch + git merge
git fetch <远程名>	        # 获取远程更新，但不自动合并，需手动检查

# 分支相关命令
git branch	                # 列出所有本地分支（当前分支前有 * 标记）
git branch -a               # 查看本地和远程所有分支
git branch <分支名>	        # 创建一个新分支
git checkout <分支名>	    # 切换到指定分支
git switch <分支名>	        # 切换到指定分支（一样，更推荐）
git merge <分支名>	        # 将指定分支合并到当前分支
git branch -d <分支名>	    # 删除已合并的本地分支
git branch -D <分支名>	    # 强制删除分支（即使未合并）
git push origin --delete <分支名>  # 删除远程分支
```
同步更新本地仓库可用（pull后的可省略）：
```bash
git fetch origin
git log origin/<分支名>
git merge <分支名>
```
或者
```bash
git pull origin <分支名>    # 拉取远程更新并自动合并到当前分支，远程名一般为origin，git pull 等同于 git fetch + git merge
```
一般上传可以用VScode的UI进行，也可以（注意替换提交说明和分支名）：
```bash
git add .
git commit -m "提交说明"
git push origin <分支名>
```
新建分支和切换如下，一般不用合并分支：
```bash
git branch <分支名>	        # 创建一个新分支
git checkout <分支名>	    # 切换到指定分支
git switch <分支名>         # 切换到指定分支（一样，更推荐）
```

# 从xml文件直接导入数据
1. 可以直接从xml文件直接读入其初始轨道参数，具体可见“<Satellite Name="Debris2" UiExpand="1">”类似字段下有如下数据：
```xml
<StartUTC>2030-11-14 08:00:00</StartUTC>
<StopUTC>2030-11-15 08:00:00</StopUTC>
<OrbEpoch>2030-11-14 08:00:00</OrbEpoch>
<StepSize>60</StepSize>
<PositionX>-1654510.742</PositionX>
<PositionY>-2551988.329</PositionY>
<PositionZ>-6763069.244</PositionZ>
<VelocityX>1342.097126</VelocityX>
<VelocityY>6593.96285</VelocityY>
<VelocityZ>-2677.689633</VelocityZ>
<GravityModel>7</GravityModel>
<MaxDegree>20</MaxDegree>
<MaxOrder>20</MaxOrder>
<UseDrag>1</UseDrag>
<UseFluxGeoFile>0</UseFluxGeoFile>
<DragCoefficient>2.2</DragCoefficient>
```
同时Satellite字段末尾会有“</Satellite>”标示结束。


# 思路
1.   先递推所有碎片位置速度序列，记录
2.   分割时间段，2小时重叠，6次机动，四段（8:00—16:00, 14:00—22:00, 20:00—4:00(+1)
和0:00(+1)—8:00(+1)）
3.   网络搜索，可以人为划定集中区间来搜索，搜索可以先粗后细，具体间隔设置多少？范可以根据碎片轨道数据定，可以覆盖其所有或者大部分轨道参数，选择优异段加入数据库（于8）
4.   可以考虑可视化所有碎片轨道参数的分布
5.   组合优化，每个组合的衔接时注意选择尽量少清除碎片的转移时间段，转移时间45min于半个周期，一般都能正常转移，还有逆行轨道，按最糟糕情况计算也可保证？（但是可能很好选出来的），组合优化按变异（随机选某一段替换），评估，交叉（随机从某一段开始叉），评估进行，覆盖局部优化（引入新基因）和全局优化（交叉，大范围组合变化），优算法的参数怎么选
6.   得到最优方案

 - 改变：如果时间太长可考虑引入束搜索，估算指标函数可根据六根数与碎片对应六根数的相对距离和数量等进行评估（未提出具体计算方案）
 - 另一种思路：网格搜索逐层收紧条件
 - J2似乎不用考虑地球自转？
# 要求
1. 判定清除碎片条件为：相对位置小于 30 km且相对速度小于 150 m/s 时
2. 母航天器采用model_rv.m递推轨道
3. 判定清除数大于5进入细搜索
4. 判定清除数大于等于8纳入优选轨道库
5. 搜索时间段分别取8:00—16:00, 14:00—22:00, 20:00—4:00(+1)和0:00(+1)—8:00(+1)，要求可指定具体时间段
6. 有疑问优先询问我进行确认

# 三母航天器组合优化

完整方法、数值轨道、18 次脉冲和验证结果见 [母航天器组合优化与机动方案](docs/母航天器组合优化与机动方案.md)。

```matlab
addpath('src');
run_mission_library_search(0); % 原有范围的三个缺失窗口粗搜诊断
run_mission_library_search(2); % 严格 >5 才细搜，>=8 入库
diagnose_reference_seeds;      % 核对参考论文的种子及角度定义
run_seeded_mission_search;     % 补充种子的局部网格细化
run_mission_combination;       % 三船联合优化、J2 两脉冲拼接、独立回放和文档导出
improve_mission_combinations;  % 同分可接受的中性变异，继续探索已有库
finalize_mission_improvements; % 高精度拼接复核补充组合，更新最终方案
crosscheck_mission_export;     % 从 CSV 重读，用 ode89 独立交叉验证
```

新结果写入 `data/mission/`，保留已有 `data/preferred_orbit_library.mat`。固定任务窗口为 UTC 08–16、14–22、20–次日04、次日00–08；其中最后两窗重叠 4 小时。最终初始状态、机动和清除清单分别为 `initial_orbits.csv`、`maneuvers.csv`、`cleared_debris.csv`。

结果为已搜索范围内的最佳可行方案，需区分本地数值复核与 ATK 官方验证；未证明全局最优。更换数据或搜索配置前，请先另存对应的 `data/mission/` 历史结果，避免入口复用旧窗口缓存。
