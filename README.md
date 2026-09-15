# ZPY15_QingZhang
第十五届全国大学生周培源力学竞赛团体赛试题方案---清障小东风队针对清华方案的复现
杨子军上传的代码

# 注意事项
1. 修改代码请先拉取远程仓库最新代码，避免冲突
2. .gitignore 文件用于忽略不需要上传的文件，如有需要请自行加入你不想上传的文件
3. 所有有单位变量注意单位，最好统一采用国际单位制，采用弧度制

# 说明
```text
  第一部分代码和第二部分混在一起
  第一部分：主要是Data_Extract_vsl,用于从ATK想定文件里提取碎片初始轨道数据并用直方图可视化统计结果。还有若干函数
  运行时需要打开ATK并加载ZPY_15.html
  第二部分：Gen_Orbit_Fragment计算24小时内碎片的位置和速度，run_all_windows是计算四个时间窗轨道基元的初筛代码，分粗筛和精筛两步（核心函数search_window.m）；refine_all_windows是用初筛的优质轨道元做种子加密M的网格，进一步精筛（核心函数refine_window.m）；命名带test的代码是测验代码。
  数据库：high_caputure_segments.mat
```