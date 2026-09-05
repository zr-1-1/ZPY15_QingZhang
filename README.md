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

# 注意事项
1. 修改代码请先拉取远程仓库最新代码，避免冲突

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