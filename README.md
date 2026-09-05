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