# ZPY15_QingZhang
第十五届全国大学生周培源力学竞赛团体赛试题方案---清障小东风队针对清华方案的复现

# 当前仓库结构
```text
ZPY15_QingZhang/
├── ZPY15.xml                         # ATK 仿真场景文件
└── src/
    ├── Data_Extract.m                # 碎片轨道数据提取与计算脚本
    ├── OE_scl_ptb.m                  # J2 摄动下的轨道根数计算
    ├── orb_elements2rv.m             # 轨道根数转位置、速度矢量
    ├── rv2coe.m                      # 位置、速度矢量转轨道根数
    └── atk_connect_matlab_dependence/ # MATLAB 与 ATK 连接所需依赖
```
`src`存放脚本，`docs`存放文档

# 常用命令
```bash
git status	# 查看工作区状态
# 同步上传相关命令
git add <文件名>	        # 将指定文件添加到暂存区
git add .	                # 添加当前目录下所有变更（新增、修改）到暂存区
git commit -m "提交说明"	# 将暂存区内容提交到本地仓库

# 分支相关命令
git branch	                # 列出所有本地分支（当前分支前有 * 标记）
git branch <分支名>	        # 创建一个新分支
git checkout <分支名>	    # 切换到指定分支
git merge <分支名>	        # 将指定分支合并到当前分支
git branch -d <分支名>	    # 删除已合并的本地分支
git branch -D <分支名>	    # 强制删除分支（即使未合并）
```