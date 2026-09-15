function benchmark_debris_grid_search_optimization
% 小规模真实数据计时，不执行完整网格，不保存/覆盖搜索结果。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'),fullfile(root,'tests'));
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
physics=struct('J_2',J_2,'R_E',R_E,'mu',mu);
timer=tic; data=load(fullfile(root,'data','x_rv.mat'),'x_rv');
fprintf('LOAD_SECONDS=%.4f\n',toc(timer));
assert(size(data.x_rv,3)>=86401,'基准需要至少 24 小时、1 秒采样的数据。');
d=pi/180;
candidates=[7000e3,.01,98*d,73*d,345*d,0; ...
    7225e3,.025,98.09*d,75.5*d,35*d,pi; ...
    7450e3,.04,98.18*d,78*d,85*d,3*pi/2];
for hours=[8,24]
    times=(0:hours*3600)';
    states=data.x_rv(:,:,1:numel(times));
    for k=1:3
        timer=tic;
        [old,oldMask,oldFirst]=reference_count_cleared_debris(candidates(k,:),states,times);
        oldSeconds=toc(timer);
        timer=tic;
        [new,newMask,newFirst]=count_cleared_debris(candidates(k,:),states,times);
        newSeconds=toc(timer);
        assert(old==new && isequal(oldMask,newMask) && isequaln(oldFirst,newFirst));
        fprintf('HOURS=%d CANDIDATE=%d OLD=%.4f NEW=%.4f CLEARED=%d\n', ...
            hours,k,oldSeconds,newSeconds,new);
    end
end
% 批量吞吐使用当前 8 小时窗口、12 个不同候选；关闭缓存以隔离并行收益。
times=(0:8*3600)'; states=data.x_rv(:,:,1:numel(times)); clear data
bounds=[7000e3,7450e3;.025,.025;98.09*d,98.09*d;75.5*d,75.5*d;35*d,35*d;0,pi];
step=[90000,1,1,1,1,pi];
opts=struct('Levels',0,'Verbose',false,'UseCache',false,'BatchSize',12);
serial=execute_debris_grid_search(states,times,bounds,step,opts, ...
    struct('EnableParallel',false),physics);
fprintf('BATCH_SERIAL_SECONDS=%.4f EVALUATIONS=%d\n',serial.Execution.SearchSeconds,serial.Evaluations);
if ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox')
    % 独立 benchmark 进程内启动 4 个进程，对比 2 和 4 个并发工作单元。
    pool=gcp('nocreate');
    if isempty(pool)
        timer=tic; pool=parpool('Processes',4);
        fprintf('POOL_START_SECONDS=%.4f\n',toc(timer));
    end
    for workers=[2,4]
        if pool.NumWorkers<workers, continue; end
        for repetition=1:2
            result=execute_debris_grid_search(states,times,bounds,step,opts, ...
                struct('EnableParallel',true,'NumWorkers',workers),physics);
            assert(isequal(result.Candidates,serial.Candidates));
            fprintf('BATCH_WORKERS=%d REP=%d SECONDS=%.4f\n', ...
                workers,repetition,result.Execution.SearchSeconds);
        end
    end
    % 数据驻留后的吞吐：同一 Constant 跨批次复用，不重复分发。
    context=struct('States',states,'Times',times,'Physics',physics);
    resident=parallel.pool.Constant(context);
    release=onCleanup(@()delete(resident)); %#ok<NASGU>
    scorer=@(oe) resident_score(oe,resident.Value);
    batch=repmat(candidates,8,1); % 24 条评分任务，不启用缓存。
    serialScorer=@(oe) count_cleared_debris(oe,states,times);
    timer=tic; expected=evaluate_grid_batch(serialScorer,batch,false);
    fprintf('RESIDENT_SERIAL_SECONDS=%.4f TASKS=%d\n',toc(timer),size(batch,1));
    for workers=[2,4]
        if pool.NumWorkers<workers, continue; end
        timer=tic;
        evaluate_grid_batch(scorer,batch(1:workers,:),true,workers);
        fprintf('RESIDENT_WARMUP_WORKERS=%d SECONDS=%.4f\n',workers,toc(timer));
        timer=tic; actual=evaluate_grid_batch(scorer,batch,true,workers);
        fprintf('RESIDENT_WORKERS=%d SECONDS=%.4f\n',workers,toc(timer));
        assert(isequal(expected,actual));
    end
end
fprintf('benchmark_debris_grid_search_optimization passed\n');
end

function score=resident_score(oe,context)
global J_2 R_E mu
J_2=context.Physics.J_2; R_E=context.Physics.R_E; mu=context.Physics.mu;
score=count_cleared_debris(oe,context.States,context.Times);
end
