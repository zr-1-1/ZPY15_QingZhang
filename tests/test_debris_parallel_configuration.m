function test_debris_parallel_configuration
% 在独立 MATLAB 测试进程运行；验证手动进程数、复用/重建和配置不落盘。
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox')
    fprintf('Parallel configuration test skipped: toolbox unavailable.\n'); return;
end
assert(isempty(gcp('nocreate')),'请在没有现有池的独立 MATLAB 进程运行此测试。');
cluster=parcluster('Processes'); originalLimit=cluster.NumWorkers;
physics=struct('J_2',0,'R_E',6378137,'mu',3.986004418e14);
global J_2 R_E mu
J_2=physics.J_2; R_E=physics.R_E; mu=physics.mu;
oe=[7200e3,.02,1,0,0,0]; states=orb_elements2rv(oe)';
bounds=[oe',oe']; step=ones(1,6);
opts=struct('Levels',0,'Verbose',false,'BatchSize',1,'ProgressIntervalSeconds',1);
settings=struct('EnableParallel',true,'NumWorkers',2,'ChunkSize',64,'WorkerThreads',1);
r=execute_debris_grid_search(states,0,bounds,step,opts,settings,physics);
assert(r.Execution.NumWorkers==2 && r.Execution.WorkerThreads==1);
assert(r.Execution.ChunkSize==64 && r.Options.BatchSize==1 && r.Options.ProgressIntervalSeconds==1);
% 本机原配置上限为 6；申请 7 个验证不再被旧上限静默截断。
settings.NumWorkers=7;
lastwarn('');
r=execute_debris_grid_search(states,0,bounds,step,opts,settings,physics);
[~,id]=lastwarn;
assert(strcmp(id,'debris:ExistingPoolLimits'));
assert(r.Execution.RequestedNumWorkers==7 && r.Execution.NumWorkers==2 && r.Execution.PoolReused);
settings.ResizeExistingPool=true;
r=execute_debris_grid_search(states,0,bounds,step,opts,settings,physics);
assert(r.Execution.NumWorkers==7 && ~r.Execution.PoolReused && r.BestScore==1);
assert(parcluster('Processes').NumWorkers==originalLimit,'不应修改磁盘上的进程池配置。');
% 调整每进程线程数也通过显式重建生效。
settings.NumWorkers=2; settings.WorkerThreads=2;
r=execute_debris_grid_search(states,0,bounds,step,opts,settings,physics);
assert(r.Execution.NumWorkers==2 && r.Execution.WorkerThreads==2);
fprintf('test_debris_parallel_configuration passed\n');
end
