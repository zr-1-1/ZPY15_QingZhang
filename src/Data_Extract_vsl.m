%% 初始化参数，与赛题说明保持一致
clc;
clear;
% 把常用参数设置为全局变量，保持一致
% J_2：J2摄动系数
% R_E：地球赤道平均半径，单位：km
% R_E_m：地球赤道平均半径，单位：m
% w_E：地球自转角速度，单位：rad/s
% mu：地球引力常数，单位：(m^3)/(s^2)
global J_2 R_E R_E_m w_E mu
% 设置参数值
C20 = -4.841653717360e-04; % 从"ATK-ZPY15专用版\AstroData\Earth\EGM96.grv"中2  0对应的完全归一化球谐系数
% Jn = -sqrt(2*n+1)*Cn0;
J_2 = -sqrt(2*2+1)*C20; % 由完全归一化球谐系数得到J2摄动系数，此处为正数
% J_2 = 1.0827e-3; % J2摄动系数，但是直接赋的值与赛题要求可能不完全一致
R_E = 6378.137; % 地球赤道平均半径，单位：km
R_E_m = 6378137; % 地球赤道平均半径，单位：m
w_E = 7.2921151467e-5; % 地球自转角速度，单位：rad/s
mu = 3.986004418e14; % 地球引力常数，单位：(m^3)/(s^2)
%% 获取所有碎片初始轨道根数
% 把connect依赖项添加到路径中
projectRoot = pwd;
% addpath(genpath(fullfile(projectRoot, 'atk_connect_matlab_dependence'))); % 把当前项目下src文件夹下的二级子目录atk_connect_matlab_dependence加入当前工作路径
addpath(genpath(projectRoot)); % 把src及其所有子目录加入当前工作路径
% 连接到atk
conID = atkOpen();
% 获取所有碎片初始轨道根数
Debris_oe = zeros(345,6);
for i = 1:345
    paramstr = sprintf('*/Satellite/Debris%d "14 Nov 2030 08:00:00.000"', i);
    Debris_d = atkConnect(conID, 'Position', paramstr);
    Debris_rv = str2double(split(Debris_d));
    Debris_oe(i,:) = rv2coe(Debris_rv(1:3),Debris_rv(4:6));
end
% 关闭atk连接C
atkClose(conID);

%% 参数可视化
titletext = ["a", "e", "i", "Omega", "omega", "M"];
figure;
tiledlayout(2,3, "Padding", "compact");
for i = 1:6
    nexttile;
    histogram(Debris_oe(:,i));
    title(titletext(i));
end