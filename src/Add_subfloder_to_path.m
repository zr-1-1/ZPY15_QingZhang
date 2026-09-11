% 把connect依赖项添加到路径中
projectRoot = pwd;
% addpath(genpath(fullfile(projectRoot, 'atk_connect_matlab_dependence'))); % 把当前项目下src文件夹下的二级子目录atk_connect_matlab_dependence加入当前工作路径
addpath(genpath(projectRoot)); % 把src及其所有子目录加入当前工作路径