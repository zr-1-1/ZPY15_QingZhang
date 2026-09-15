%% 独立入口：使用已保存的碎片状态执行粗搜—细搜，无需连接 ATK。
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'src'));
global J_2 R_E mu
J_2 = -sqrt(5)*(-4.841653717360e-04);
R_E = 6378137;
mu = 3.986004418e14;
d2r = pi/180;

%% 性能配置（i5-13600KF / 14 核 20 线程 / 32 GB 内存的保守起点）
% 所有性能参数集中在这里；工作进程数不等于 CPU 逻辑线程数。
performance.EnableParallel = true;
performance.NumWorkers = 6;           % 可手动改为 4、8、10 等；新建池按此值申请。
performance.WorkerThreads = 1;        % 每个 MATLAB 工作进程的计算线程数，避免过度争用。
performance.ChunkSize = 512;          % 清除判定时间块，单位：帧；不改变采样间隔。
performance.BatchSize = 128;          % 每批候选轨道数。
performance.UseCache = true;
performance.ProgressIntervalSeconds = 10; % 完成批次后检查是否需要更新进度。
% performance.ResizeExistingPool = false; % true：允许重建空闲池以匹配进程数/线程数。
performance.ResizeExistingPool = true; % true：允许重建空闲池以匹配进程数/线程数。
% false 时复用已有池，明确报告请求/实际配置；不会中断已有后台任务。

%% 可调参数（a 为 m，角度为 rad；根数范围对应 searchStartUTC 时刻）
% 以下为 UTC，不是北京时间。数据历元必须与生成 x_rv 时一致。
stateDataStartUTC = datetime(2030,11,14,8,0,0,'TimeZone','UTC');
searchStartUTC = datetime(2030,11,14,8,0,0,'TimeZone','UTC');
searchEndUTC = datetime(2030,11,14,16,0,0,'TimeZone','UTC');
stateStepSeconds = 1; % 必须与生成 x_rv 时的采样间隔一致。
searchBounds = [7000e3,7450e3; .01,.04; 98*d2r,98.18*d2r; ...
    73*d2r,78*d2r; 345*d2r,445*d2r; 0,2*pi];
coarseStep = [150e3,.015,.09*d2r,2.5*d2r,50*d2r,90*d2r];
target_min = 3; % 清除数严格大于 target_min 才进入细搜范围。
target_num = 8; % 清除数 >= target_num 的所有已评估轨道存入优选轨道库。
searchOptions = struct('TopK',3,'Levels',2,'RefineFactor',2, ...
    'MaxEvaluations',200000,'Verbose',true,'RefineThreshold',target_min, ...
    'ArchiveThreshold',target_num,'BatchSize',performance.BatchSize, ...
    'UseCache',performance.UseCache,'ProgressIntervalSeconds',performance.ProgressIntervalSeconds);
% 有并行工具箱时使用进程池；不可用时自动回退串行，不改变搜索规则。
executionOptions = struct('EnableParallel',performance.EnableParallel, ...
    'NumWorkers',performance.NumWorkers,'WorkerThreads',performance.WorkerThreads, ...
    'ChunkSize',performance.ChunkSize,'ResizeExistingPool',performance.ResizeExistingPool);

%% 读取数据；保留每个时间采样点，粗搜与细搜的时间分辨率一致。
stateFile = fullfile(projectRoot,'data','x_rv.mat');
assert(isfile(stateFile),'请先生成 data/x_rv.mat，再运行本脚本。');
loadedStates = load(stateFile,'x_rv');
assert(isfield(loadedStates,'x_rv'),'MAT 文件中缺少 x_rv。');
x_rv = loadedStates.x_rv;
clear loadedStates
validateattributes(x_rv,{'double','single'},{'real','nonempty','finite'});
assert(size(x_rv,1)==6 && ndims(x_rv)<=3,'x_rv 必须为 6 x N x T。');
validateattributes(stateStepSeconds,{'numeric'},{'scalar','positive','finite'});
[searchFrameIndices,time] = debris_search_time_window(stateDataStartUTC, ...
    size(x_rv,3),stateStepSeconds,searchStartUTC,searchEndUTC);
x_rv = x_rv(:,:,searchFrameIndices);
fprintf('搜索时间（UTC）：%s 至 %s，共 %d 帧。\n', ...
    char(searchStartUTC),char(searchEndUTC),numel(time));
physics = struct('J_2',J_2,'R_E',R_E,'mu',mu);
searchResult = execute_debris_grid_search(x_rv,time,searchBounds,coarseStep, ...
    searchOptions,executionOptions,physics);
[~,searchResult.ClearedMask,searchResult.FirstClearTime] = ...
    count_cleared_debris(searchResult.BestOE,x_rv,time,executionOptions.ChunkSize);
searchResult.ClearedIDs = find(searchResult.ClearedMask);
searchResult.EligibleForRefinement = searchResult.BestScore>searchOptions.RefineThreshold;
searchResult.TargetMin = target_min;
searchResult.TargetNum = target_num;
searchResult.MeetsTarget = searchResult.BestScore>=target_num;
searchResult.StateFile = stateFile;
searchResult.StateStepSeconds = stateStepSeconds;
searchResult.StateDataStartUTC = stateDataStartUTC;
searchResult.SearchStartUTC = searchStartUTC;
searchResult.SearchEndUTC = searchEndUTC;
searchResult.OrbitEpochUTC = searchStartUTC;
searchResult.FirstClearUTC = searchStartUTC+seconds(searchResult.FirstClearTime);
searchResult.DistanceThreshold = 30000;
searchResult.SpeedThreshold = 150;
searchResult.PropagationModel = 'ode45 + model_rv';
save(fullfile(projectRoot,'data','grid_search_result.mat'),'searchResult');
preferredOrbitLibrary = struct('Orbits',searchResult.PreferredOrbits, ...
    'Columns',{{'a_m','e','i_rad','Omega_rad','omega_rad','M_rad','ClearedCount'}}, ...
    'TargetNum',target_num,'StateFile',stateFile, ...
    'StateStepSeconds',stateStepSeconds,'StateDataStartUTC',stateDataStartUTC, ...
    'SearchStartUTC',searchStartUTC,'SearchEndUTC',searchEndUTC, ...
    'OrbitEpochUTC',searchStartUTC, ...
    'PropagationModel',searchResult.PropagationModel, ...
    'DistanceThreshold',30000,'SpeedThreshold',150);
save(fullfile(projectRoot,'data','preferred_orbit_library.mat'),'preferredOrbitLibrary');
fprintf('优选轨道库：%d 条清除数 >= %d 的轨道。\n', ...
    size(preferredOrbitLibrary.Orbits,1),target_num);
fprintf('最高采样清除数：%d；目标 %d；结果已保存。\n',searchResult.BestScore,target_num);
disp('最优根数 [a,e,i,Omega,omega,M]（m、rad）：');
disp(searchResult.BestOE);
disp('已清除碎片编号（x_rv 的列编号）：');
disp(searchResult.ClearedIDs);
