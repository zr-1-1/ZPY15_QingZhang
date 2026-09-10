%% 初始化参数，与赛题说明保持一致
clc;
clear;
d2r = pi/180; % 角度转弧度
r2d = 180/pi; % 弧度转角度
% 把常用参数设置为全局变量，保持一致
% J_2：J2摄动系数
% R_E：地球赤道平均半径，单位：m
% w_E：地球自转角速度，单位：rad/s
% mu：地球引力常数，单位：(m^3)/(s^2)
% 轨道根数顺序：[a(m), e, i(rad), Omega(rad), omega(rad), M(rad)]
global J_2 R_E w_E mu
% 设置参数值
C20 = -4.841653717360e-04; % 从"ATK-ZPY15专用版\AstroData\Earth\EGM96.grv"中2  0对应的完全归一化球谐系数
% Jn = -sqrt(2*n+1)*Cn0;
J_2 = -sqrt(2*2+1)*C20; % 由完全归一化球谐系数得到J2摄动系数，此处为正数
% J_2 = 1.0827e-3; % J2摄动系数，但是直接赋的值与赛题要求可能不完全一致
R_E = 6378137; % 地球赤道平均半径，单位：m
w_E = 7.2921151467e-5; % 地球自转角速度，单位：rad/s
mu = 3.986004418e14; % 地球引力常数，单位：(m^3)/(s^2)
projectRoot = pwd;
%% 获取所有碎片初始轨道根数
% 把connect依赖项添加到路径中
projectRoot = pwd;
% addpath(genpath(fullfile(projectRoot, "src", 'atk_connect_matlab_dependence'))); % 把当前项目下src文件夹下的二级子目录atk_connect_matlab_dependence加入当前工作路径
addpath(genpath(fullfile(projectRoot, "src"))); % 把src及其所有子目录加入当前工作路径
% 连接到atk
conID = atkOpen();
% 获取所有碎片初始轨道根数
Debris_oe = zeros(345,6);
for i = 1:345
    paramstr = sprintf('*/Satellite/Debris%d "14 Nov 2030 08:00:00.000"', i);
    Debris_d = atkConnect(conID, 'Position', paramstr);
    Debris_rv = str2double(split(Debris_d)); % 接口输入约定：位置 m、速度 m/s，需与 ATK 实际输出单位一致
    Debris_oe(i,:) = rv2coe(Debris_rv(1:3),Debris_rv(4:6));
end
% 保存所有碎片初始轨道根数对应变量Debris_oe为静态文件，之后不用再从ATK导入
save(fullfile(projectRoot,'data','data.mat'), 'Debris_oe'); % 该变量不带对应时间信息，但是所有碎片初始轨道根数历元都为最初时刻2030-11-14 08:00:00.000(UTC)
% 通过load('data.mat', 'Debris_oe');加载该变量到matlab工作区
% 关闭atk连接C
atkClose(conID);

%% 按照一定步长生成24小时内碎片每个时刻的轨道根数
load(fullfile(projectRoot,'data','data.mat'), 'Debris_oe');
%% 直接读取从xml文件得到的数据
debris_orbits = readtable(fullfile(projectRoot,'data','debris_orbits.csv'));
num = height(debris_orbits);
% 可用debris_orbits(1,1)读取导入csv对应行列的数据
% 第5列为历元，第7~12列为位置和速度三轴分量(单位为m和m/s)（ICRF地心惯性系下的三轴分量）
% debris_orbits.PositionX(3);debris_orbits.VelocityY(6);分别读取第三个碎片X轴位置和第六个碎片Y轴速度分量
Debris_rv = zeros(6, num);
Debris_oe = zeros(6, num);
for i = 1:num
    Debris_rv(:,i) = [
    debris_orbits.PositionX(i);
    debris_orbits.PositionY(i);
    debris_orbits.PositionZ(i);
    debris_orbits.VelocityX(i);
    debris_orbits.VelocityY(i);
    debris_orbits.VelocityZ(i);
    ];
    temp = rv2orb_elements(Debris_rv(1:3,i),Debris_rv(4:6,i), mu);
    Debris_oe(:,i) = temp;
end
Debris_oe = Debris_oe';

% load(fullfile(projectRoot,'data','debris_initial_state.mat'), 'debris_initial_state');
%% 轨道根数可视化
% 只转换显示副本，Debris_oe 仍保持 m、rad，供后续计算使用。
assert(size(Debris_oe,2) == 6 && ~isempty(Debris_oe) && ...
    isreal(Debris_oe) && all(isfinite(Debris_oe),'all'), ...
    '轨道根数必须是非空、有限实数的 N×6 矩阵。');
oeNames = ["a 半长轴"; "e 偏心率"; "i 倾角"; ...
    "Omega 升交点赤经"; "omega 近地点幅角"; "M 平近点角"];
oeUnits = ["km"; "无量纲"; "deg"; "deg"; "deg"; "deg"];
oeDisplayScale = [1e-3, 1, 180/pi, 180/pi, 180/pi, 180/pi];
oeDisplay = Debris_oe .* oeDisplayScale;
% 周期角沿最大空白间隔切开，得到包含全部样本的最短连续圆弧。
% 跨越 0° 时允许上界超过 360°，例如 [350,370] 表示 [350,360)∪[0,10]。
for oeIdx = 4:6
    oeAngles = sort(mod(oeDisplay(:,oeIdx),360));
    [~,oeGapIdx] = max(diff([oeAngles; oeAngles(1)+360]));
    oeArcStart = oeAngles(mod(oeGapIdx,numel(oeAngles))+1);
    oeWrapped = mod(oeDisplay(:,oeIdx),360);
    oeWrapped(oeWrapped < oeArcStart) = oeWrapped(oeWrapped < oeArcStart)+360;
    oeDisplay(:,oeIdx) = oeWrapped;
end
oeSampleMin = min(oeDisplay,[],1);
oeSampleMax = max(oeDisplay,[],1);
% 可调整：只控制边界向外取整，不是网格步长；单位与图中一致。
oeBoundaryRound = [100, 0.01, 0.5, 5, 5, 5];
oeGridLower = floor(oeSampleMin./oeBoundaryRound).*oeBoundaryRound;
oeGridUpper = ceil(oeSampleMax./oeBoundaryRound).*oeBoundaryRound;
% 单值样本也提供非零区间；限制非周期根数的物理范围。
oeSameBound = oeGridLower == oeGridUpper;
oeGridUpper(oeSameBound) = oeGridUpper(oeSameBound)+oeBoundaryRound(oeSameBound);
oeGridLower(1:3) = max(oeGridLower(1:3),0);
oeGridUpper(2:3) = min(oeGridUpper(2:3),[1,180]);
for oeIdx = 4:6
    if oeGridUpper(oeIdx)-oeGridLower(oeIdx) >= 360
        oeGridLower(oeIdx) = 0;
        oeGridUpper(oeIdx) = 360;
        oeDisplay(:,oeIdx) = mod(oeDisplay(:,oeIdx),360);
        oeSampleMin(oeIdx) = min(oeDisplay(:,oeIdx));
        oeSampleMax(oeIdx) = max(oeDisplay(:,oeIdx));
    end
end
% 六行分别对应 a,e,i,Omega,omega,M，两列为下/上界，单位恢复为 m、rad。
oeGridBoundsSI = [oeGridLower(:),oeGridUpper(:)]./oeDisplayScale(:);
oeRangeTable = table(oeNames,oeUnits,oeSampleMin(:),oeSampleMax(:), ...
    oeGridLower(:),oeGridUpper(:),(oeGridUpper-oeGridLower)', ...
    'VariableNames',{'Element','Unit','SampleMin','SampleMax', ...
    'GridLower','GridUpper','GridWidth'});
disp('初始样本覆盖范围（向外取整；不是最优解范围或网格步长）：');
disp(oeRangeTable);
fprintf(['计算边界见 oeGridBoundsSI（m、rad）。周期角采样后用 mod(angle,2*pi)，' ...
    '整圈网格不要重复包含首尾点。\n']);

% 每次运行本节，将当前范围以 UTF-8 文本保存到 data 文件夹。
oeRangeFile = fullfile(projectRoot,'data','初始样本覆盖范围.txt');
if ~isfolder(fileparts(oeRangeFile))
    mkdir(fileparts(oeRangeFile));
end
oeRangeText = ["初始样本覆盖范围"; ...
    "生成时间：" + string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')); ...
    "样本数量：" + string(size(Debris_oe,1)); ...
    "说明：范围仅覆盖当前历元样本；搜索参考边界向外取整，不代表最优解范围。"; ...
    "周期角采用连续展开区间，超过 360 deg 与减去 360 deg 等价。"; ...
    "边界取整单位不是网格步长；整圈角度网格不要重复包含首尾点。"; ""];
for oeIdx = 1:6
    oeRangeText(end+1,1) = sprintf('%s（%s）',oeNames(oeIdx),oeUnits(oeIdx));
    oeRangeText(end+1,1) = sprintf('  样本范围：[%.15g, %.15g]', ...
        oeSampleMin(oeIdx),oeSampleMax(oeIdx));
    oeRangeText(end+1,1) = sprintf('  搜索参考：[%.15g, %.15g]；宽度：%.15g；边界取整单位：%.15g', ...
        oeGridLower(oeIdx),oeGridUpper(oeIdx), ...
        oeGridUpper(oeIdx)-oeGridLower(oeIdx),oeBoundaryRound(oeIdx));
    oeRangeText(end+1,1) = "";
end
oeRangeText(end+1,1) = "计算用搜索参考边界（a: m；e: 无量纲；角度: rad）：";
for oeIdx = 1:6
    oeRangeText(end+1,1) = sprintf('%s：[%.15g, %.15g]', ...
        oeNames(oeIdx),oeGridBoundsSI(oeIdx,1),oeGridBoundsSI(oeIdx,2));
end
writelines(oeRangeText,oeRangeFile,'Encoding','UTF-8');
fprintf('初始样本覆盖范围已保存：%s\n',oeRangeFile);

oeFigure = figure('Name','碎片轨道根数与搜索范围','Color','w', ...
    'Position',[80,80,1380,820]);
oeLayout = tiledlayout(oeFigure,2,3,'Padding','compact','TileSpacing','loose');
title(oeLayout,sprintf('%d 个碎片 | 初始轨道根数与搜索范围',size(Debris_oe,1)));
subtitle(oeLayout,{'蓝色：样本分布；橙色虚线：向外取整边界（不是搜索步长）', ...
    '周期角按最短覆盖圆弧展开，超过 360° 与减去 360° 等价；范围仅覆盖当前历元样本'});
for oeIdx = 1:6
    oeAxes = nexttile(oeLayout);
    histogram(oeAxes,oeDisplay(:,oeIdx),'BinMethod','fd', ...
        'BinLimits',[oeGridLower(oeIdx),oeGridUpper(oeIdx)], ...
        'FaceColor',[0.18,0.46,0.70],'FaceAlpha',0.85,'EdgeColor','w');
    hold(oeAxes,'on');
    xline(oeAxes,oeGridLower(oeIdx),'--','Color',[0.85,0.35,0.10],'LineWidth',1.5);
    xline(oeAxes,oeGridUpper(oeIdx),'--','Color',[0.85,0.35,0.10],'LineWidth',1.5);
    oePad = max(oeGridUpper(oeIdx)-oeGridLower(oeIdx),oeBoundaryRound(oeIdx))*0.08;
    xlim(oeAxes,[oeGridLower(oeIdx)-oePad,oeGridUpper(oeIdx)+oePad]);
    title(oeAxes,{char(oeNames(oeIdx)), ...
        sprintf('样本 [%.6g, %.6g]',oeSampleMin(oeIdx),oeSampleMax(oeIdx)), ...
        sprintf('搜索参考 [%.6g, %.6g]',oeGridLower(oeIdx),oeGridUpper(oeIdx))}, ...
        'Interpreter','none','FontSize',11);
    xlabel(oeAxes,oeUnits(oeIdx));
    ylabel(oeAxes,'碎片数');
    grid(oeAxes,'on');
    oeAxes.GridAlpha = 0.15;
    oeAxes.FontSize = 10;
    oeAxes.XAxis.Exponent = 0;
end

%% 网格搜索，先粗后细
step = [10*1000 0.001 0.01*d2r 0.01*d2r 1*d2r 1*d2r]'; % 粗搜索间隔，单位为m和弧度
range = [7000*1000 7443.24*1000
         0.00001 0.0308333
         98*d2r 98.1705*d2r
         73.8298*d2r 77.6596*d2r
         345*d2r 445*d2r
         0 360*d2r
];
target_min = 4;
target_num = 8;


%% 
t_step = 1; % 步长，单位：秒
outercell = cell(1,345);
for i = 1:345
    innercell = cell(1,20000);
    E_0 = Debris_oe(i, :);
    idx = 0;

    for t = 0: t_step: 60*60*24
        E_t = OE_scl_ptb(E_0, t);
        idx = idx+1;
        innercell{idx} = E_t;
    end
    D_E_i = cat(2,innercell{:});
    outercell{i} = D_E_i;
end
D_Eall = cat(3,outercell{:});

%% 子函数
function [orbital_elements] = rv2orb_elements(r, v, mu)
% r为初始位置矢量(列向量)，单位：m
% v为初始速度矢量(列向量)，单位：m/s
% miu为地球引力常数，单位：m^3/s^2
% 暂时未做圆轨道和赤道轨道的指定和兼容
% 对当前碎片轨道根数转换影响不大

% 强制转为列向量
r = r(:); v = v(:);
h = cross(r, v); % 计算轨道角动量矢量，单位：m^2/s
h_mo = sqrt(h' * h); % 计算轨道角动量大小，单位：m^2/s
cosi = h(3)/h_mo; % 计算轨道倾角余弦值
omega = atan2(h(1), -h(2)); % 计算升交点赤经，单位：rad
i = acos(cosi); % 计算轨道倾角，单位：rad
% i_deg = i / pi * 180; % 计算轨道倾角，单位：°
% omega_deg = omega / pi * 180; % 计算升交点赤经，单位：°

e_vector = (cross(v, h)/mu) - (r/norm(r)); % 计算偏心率矢量
e = norm(e_vector); % 计算轨道偏心率

sinw = e_vector(3)/(e*sin(i)); % 计算近地点幅角正弦值
cosw = (e_vector(1)*cos(omega)+e_vector(2)*sin(omega))/(e); % 计算近地点幅角余弦值
w = atan2(sinw, cosw); % 计算近地点幅角，单位：rad
% w_deg = w / pi * 180; % 计算近地点幅角，单位：°
a = h_mo^2/(mu*(1-e^2)); % 计算轨道半长轴，单位：m
sin_theta_w = r(3)/(norm(r)*sin(i));  % 计算真近点角正弦值
cos_theta_w = (r(1)*cos(omega)+r(2)*sin(omega))/(norm(r)); % 计算真近点角余弦值
theta = atan2(sin_theta_w, cos_theta_w) - w; % 计算真近点角，单位：rad
% 或者用位置和偏心率矢量的有向夹角求真近点角：
% cos_theta = dot(e_vector,r)/(e*norm(r));
% sin_theta = dot(cross(e_vector,r),h)/(e*norm(r)*norm(h));
% theta = atan2(sin_theta,cos_theta);
cos_phi = (e*a+norm(r)*cos(theta))/a; % 偏近点角余弦值
sin_phi = sqrt(1-e^2)*sin(theta)/(1+e*cos(theta)); % 偏近点角正弦值
phi = atan2(sin_phi, cos_phi); % 计算偏近点角，单位：rad
M = phi - e*sin(phi); % 计算平近点角，单位：rad
M = mod(M,2*pi); % 约束平近点角范围到[0, 2π)
% theta_deg = theta / pi * 180; % 计算真近点角，单位：°
% r'*v/(norm(r)*norm(v)); % 计算径向速度方向余弦值
orbital_elements = [a;e;i;omega;w;M];
end

function [r,v] = orbit6toRV(orbit_elements, miu)
% 轨道六根数转地心惯性系下的速度位置矢量的函数
a = orbit_elements(1);       % 半长轴(m)
e = orbit_elements(2);       % 偏心率
i = orbit_elements(3);       % 轨道倾角(rad)
omega = orbit_elements(4);   % 升交点赤经(rad)
w = orbit_elements(5);       % 近地点幅角(rad)
M = orbit_elements(6);   % 平近点角(rad)
% theta = orbit_elements(6);   % 真近点角(rad)
E = solve_Kepler(M, e); % 计算偏近点角(rad)
theta = E2theta(E,e);   % 计算真近点角(rad)
r_dist = a*(1-e^2)/(1+e*cos(theta)); % 地心距离
r_o = [r_dist*cos(theta);r_dist*sin(theta);0]; % 轨道系下地心矢量
C1 = EulerZ2C(-w);
C2 = EulerX2C(-i);
C3 = EulerZ2C(-omega);
r_ei = C3*C2*C1*r_o; % 地心惯性系下的地心矢量(m)
r = r_ei;                % 位置矢量(m)

% 直接按照课程ppt里的公式计算
u = w + theta; % 纬度幅角(rad)
v0 = [-cos(omega)*(sin(u)+e*sin(w))-sin(omega)*(cos(u)+e*cos(w))*cos(i)
      -sin(omega)*(sin(u)+e*sin(w))+cos(omega)*(cos(u)+e*cos(w))*cos(i)
      (cos(u)+e*cos(w))*sin(i)]; % 速度在地心惯性系下的方向矢量（不是单位矢量）
v = sqrt(miu/(a*(1-e^2))).*v0;   % 在地心惯性系下的速度矢量(m/s)
end

function [x] = sat_sim(r0,v0,Tr,ts)
% 母航天器仿真函数，Tr为仿真总时间，ts为仿真步长(0~Tr)，输入可整除的数
% 输出x为列向量组，前三排为r，后三排为v
% 强制转为列向量
r0 = r0(:); v0 = v0(:);
n=round(Tr/ts);          % 仿真步数
% time = (0 : ts : Tr)';  % 仿真时间向量，包含 t = 0 到 t = Tr，共 n+1 个点
n0 = 6;             % 状态变量个数
x = zeros(n0,n+1);   % 状态向量
x(:,1) = [r0;v0];
    for i=1:1:n
        tSpan=[0 ts];
        [tt,xx] = ode45(@model_rv,tSpan,x(:,i),[]);  % 求解微分方程
        x(:,i+1) = xx(end,:)'; % 记录状态变量
    end
end

function E = solve_Kepler(M, e)
% 用牛顿法求解开普勒方程的函数
% 输入平近点角M(rad)和偏心率e，输出偏近点角E(rad)
d2r = pi/180; % 角度转弧度
% 使用牛顿-辛普森方法求解
E = M;
% max_times = 1000;
% i = 0;
while(1)
    E_last = E;
    E = E - (E-e*sin(E)-M)/(1-e*cos(E));
    if abs(E - E_last) < 0.001*d2r
        break;
    end
end
end

function theta = E2theta(E,e)
% 偏近点角(rad)转真近点角(rad),0<e<1
if E == pi
    theta = pi;
else
    theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    % 约束范围到[0,2*pi]
    % if theta<0
    % theta = theta + 2*pi;
    % end
    theta = mod(theta,2*pi); % 约束范围到[0, 2π)
end
end

function E = theta2E(theta,e)
% 真近点角(rad)转偏近点角(rad),0<e<1
if theta == pi
    E = pi;
else
    E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
    % 约束范围到[0,2*pi]
    % if E<0
    % E = E + 2*pi;
    % endE
    E = mod(E,2*pi); % 约束范围到[0, 2π)
end
end

function Cz = EulerZ2C(EulerZ)
% 由z轴转角得到对应的基本方向余弦矩阵
Cz = [cos(EulerZ) sin(EulerZ) 0
    -sin(EulerZ) cos(EulerZ) 0
    0 0 1];
end

function Cy = EulerY2C(EulerY)
% 由y轴转角得到对应的基本方向余弦矩阵
Cy = [cos(EulerY) 0 -sin(EulerY)
    0 1 0
    sin(EulerY) 0 cos(EulerY)];
end

function Cx = EulerX2C(EulerX)
% 由x轴转角得到对应的基本方向余弦矩阵
Cx = [1 0 0
    0 cos(EulerX) sin(EulerX)
    0 -sin(EulerX) cos(EulerX)];
end