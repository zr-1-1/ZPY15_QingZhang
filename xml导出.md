
## 导出为 CSV（无需启动 ATK）

在仓库根目录运行（需要 Python 3，无需安装第三方库）：

```bash
python src/export_debris_csv.py
# 可选：指定输入文件和输出位置
python src/export_debris_csv.py ZPY15.xml -o data/debris_orbits.csv
```

默认在仓库根目录生成 `debris_orbits.csv`，再次运行会覆盖该文件。默认输入、输出路径相对于脚本所在仓库定位，从其他工作目录调用也可用；手动指定的相对路径相对于当前工作目录。XML 按文件声明的 GB2312 编码读取；其他编码的文件可用 `--encoding utf-8` 等显式指定。

当前场景导出 345 行，每个碎片一行，按 `DebrisID` 数值排序（Debris1、Debris2、…、Debris345）。仅读取 `Satellite/Orbit` 的直接字段，避免混入第三天体引力配置中的同名字段。列含义如下：

| 列 | 含义与单位 |
| --- | --- |
| `DebrisID`, `Name` | 碎片编号、XML 中的名称 |
| `StartUTC`, `StopUTC`, `OrbEpoch` | 仿真起止时间、初始状态历元，均为 UTC |
| `StepSize` | XML 中的轨道步长，秒 |
| `PositionX`, `PositionY`, `PositionZ` | 历元时刻的位置，米 |
| `VelocityX`, `VelocityY`, `VelocityZ` | 历元时刻的速度，米/秒 |
| `GravityModel`, `MaxDegree`, `MaxOrder` | 原始引力模型编号及阶次配置 |
| `UseDrag`, `UseFluxGeoFile`, `DragCoefficient` | 原始阻力开关、通量文件开关和阻力系数 |

数值保留 XML 中的文本精度，不进行单位或坐标系转换。CSV 保存的是**初始笛卡尔状态及上述配置**，不是六个轨道根数，也不是按步长传播后的轨迹；它不包含完整的 ATK 摄动模型配置。

MATLAB 在仓库根目录读取：

```matlab
T = readtable('debris_orbits.csv');
debris_id = T.DebrisID;
r0 = T{:, {'PositionX', 'PositionY', 'PositionZ'}}; % N×3，m
v0 = T{:, {'VelocityX', 'VelocityY', 'VelocityZ'}}; % N×3，m/s
epoch = datetime(T.OrbEpoch, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
                 'TimeZone', 'UTC');
save('debris_initial_state.mat', 'debris_id', 'r0', 'v0', 'epoch');
```

以后直接 `load('debris_initial_state.mat')` 即可。如果后续计算使用 SI 制引力常数 `mu = 3.986004418e14`，应直接使用这里的米和米/秒数据，不要再除以 1000。

## 更简单的方案：MATLAB 直接读 XML

如果只在 MATLAB 中使用，可以跳过 Python 和 CSV，用 `xmlread` 提取初始状态，再保存为 `.mat`。以下示例在仓库根目录运行，按编号排序：

```matlab
root = pwd;
doc = xmlread(fullfile(root, 'ZPY15.xml'));
satellites = doc.getDocumentElement().getElementsByTagName('Satellite');
fields = {'PositionX', 'PositionY', 'PositionZ', ...
          'VelocityX', 'VelocityY', 'VelocityZ'};
state = zeros(0, 6);
debris_id = zeros(0, 1);
epoch_text = cell(0, 1);
for k = 0:satellites.getLength()-1
    sat = satellites.item(k);
    token = regexp(char(sat.getAttribute('Name')), '^Debris(\d+)$', 'tokens', 'once');
    if isempty(token), continue; end
    orbit = sat.getElementsByTagName('Orbit').item(0);
    row = numel(debris_id) + 1;
    debris_id(row, 1) = str2double(token{1});
    epoch_text{row, 1} = char(orbit.getElementsByTagName('OrbEpoch').item(0).getTextContent());
    for j = 1:6
        node = orbit.getElementsByTagName(fields{j}).item(0);
        state(row, j) = str2double(char(node.getTextContent()));
    end
end
[debris_id, order] = sort(debris_id);
r0 = state(order, 1:3); % m
v0 = state(order, 4:6); % m/s
epoch = datetime(epoch_text(order), 'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
                 'TimeZone', 'UTC');
save(fullfile(root,'data','debris_initial_state.mat'), 'debris_id', 'r0', 'v0', 'epoch');
```
