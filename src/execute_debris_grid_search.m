function result = execute_debris_grid_search(debrisStates, times, bounds, step, searchOptions, executionOptions, physics)
% 配置串行/进程池评分。每个工作进程只接收一份时间窗口内的碎片数据。
% executionOptions: EnableParallel(true), NumWorkers(6), WorkerThreads(1),
% ChunkSize(512), ResizeExistingPool(false)。仅修改本次 cluster 对象，不保存配置文件。
% physics 含 J_2、R_E、mu；每次评分显式初始化工作进程的全局参数。
defaults=struct('EnableParallel',true,'NumWorkers',6,'ChunkSize',512, ...
    'WorkerThreads',1,'ResizeExistingPool',false);
names=fieldnames(executionOptions);
for k=1:numel(names)
    assert(isfield(defaults,names{k}),'未知执行选项：%s',names{k});
    defaults.(names{k})=executionOptions.(names{k});
end
validateattributes(defaults.EnableParallel,{'logical'},{'scalar'});
validateattributes(defaults.NumWorkers,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(defaults.ChunkSize,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(defaults.WorkerThreads,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(defaults.ResizeExistingPool,{'logical'},{'scalar'});
validateattributes(physics.J_2,{'numeric'},{'scalar','real','finite'});
validateattributes(physics.R_E,{'numeric'},{'scalar','positive','finite'});
validateattributes(physics.mu,{'numeric'},{'scalar','positive','finite'});
execution=struct('ParallelRequested',defaults.EnableParallel,'ParallelUsed',false, ...
    'NumWorkers',0,'FallbackReason','','ChunkSize',defaults.ChunkSize, ...
    'PoolSetupSeconds',0,'RequestedNumWorkers',defaults.NumWorkers, ...
    'RequestedWorkerThreads',defaults.WorkerThreads,'WorkerThreads',0, ...
    'PoolNumWorkers',0,'PoolReused',false);
poolTimer=tic;
if defaults.EnableParallel
    try
        assert(~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox'), ...
            'Parallel Computing Toolbox 不可用。');
        pool=gcp('nocreate');
        if ~isempty(pool) && isa(pool,'parallel.ProcessPool') && ...
                defaults.ResizeExistingPool && ...
                (pool.NumWorkers~=defaults.NumWorkers || pool.Cluster.NumThreads~=defaults.WorkerThreads)
            % 上次 Constant 释放可能留下短暂的异步清理任务，先给它们结束的机会。
            idleTimer=tic;
            while (~isempty(pool.FevalQueue.RunningFutures) || ~isempty(pool.FevalQueue.QueuedFutures)) ...
                    && toc(idleTimer)<2
                pause(0.05);
            end
            assert(isempty(pool.FevalQueue.RunningFutures) && isempty(pool.FevalQueue.QueuedFutures), ...
                '现有池有后台任务，不能重建；保留现有池并回退串行。');
            delete(pool);
            pool=[];
        end
        execution.PoolReused=~isempty(pool);
        if isempty(pool)
            cluster=parcluster('Processes');
            % 原配置上限可能是 6：只调整内存中的对象，让用户指定的 8/10 等生效。
            cluster.NumWorkers=max(cluster.NumWorkers,defaults.NumWorkers);
            cluster.NumThreads=defaults.WorkerThreads;
            pool=parpool(cluster,defaults.NumWorkers);
        end
        assert(isa(pool,'parallel.ProcessPool'), ...
            '现有池不是进程池；保留现有池并使用串行评分。');
        execution.ParallelUsed=true;
        execution.NumWorkers=min(defaults.NumWorkers,pool.NumWorkers);
        execution.PoolNumWorkers=pool.NumWorkers;
        execution.WorkerThreads=pool.Cluster.NumThreads;
        if execution.NumWorkers~=defaults.NumWorkers || execution.WorkerThreads~=defaults.WorkerThreads
            warning('debris:ExistingPoolLimits', ...
                ['请求 %d 个进程、每进程 %d 线程；现有池实际使用 %d 个进程、每进程 %d 线程。' ...
                '可设置 ResizeExistingPool=true，在无后台任务时重建池。'], ...
                defaults.NumWorkers,defaults.WorkerThreads,execution.NumWorkers,execution.WorkerThreads);
        end
    catch failure
        execution.FallbackReason=failure.message;
        warning('debris:ParallelUnavailable','并行初始化失败，改用串行评分：%s',failure.message);
    end
end
execution.PoolSetupSeconds=toc(poolTimer);
context=struct('States',debrisStates,'Times',times,'Physics',physics, ...
    'ChunkSize',defaults.ChunkSize);
if execution.ParallelUsed
    data=parallel.pool.Constant(context);
    releaseData=onCleanup(@()delete(data)); %#ok<NASGU>
    scoreFcn=@(oe) score_debris_context(oe,data.Value);
else
    scoreFcn=@(oe) score_debris_context(oe,context);
end
searchOptions.UseParallel=execution.ParallelUsed;
searchOptions.MaxWorkers=max(1,execution.NumWorkers);
if ~isfield(searchOptions,'Verbose') || searchOptions.Verbose
    fprintf('评分方式：并行=%d，请求/实际进程=%d/%d，每进程线程=%d，时间块=%d 帧。\n', ...
        execution.ParallelUsed,defaults.NumWorkers,execution.NumWorkers, ...
        execution.WorkerThreads,defaults.ChunkSize);
end
searchTimer=tic;
result=coarse_to_fine_grid_search(scoreFcn,bounds,step,searchOptions);
execution.SearchSeconds=toc(searchTimer);
result.Execution=execution;
end

function score = score_debris_context(oe,context)
global J_2 R_E mu
J_2=context.Physics.J_2;
R_E=context.Physics.R_E;
mu=context.Physics.mu;
score=count_cleared_debris(oe,context.States,context.Times,context.ChunkSize);
end
