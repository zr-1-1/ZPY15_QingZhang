function test_grid_search_optimization
% 优化前后结果对比；使用真实动力学和合成碎片，不写搜索结果文件。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'),fullfile(root,'tests'));
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
physics=struct('J_2',J_2,'R_E',R_E,'mu',mu);
oe=[7200e3,.025,1.71,1.3,.4,.6];
times=(0:1024)';
sol=ode45(@model_rv,[0,times(end)],orb_elements2rv(oe)', ...
    odeset('RelTol',1e-9,'AbsTol',[1e-3;1e-3;1e-3;1e-6;1e-6;1e-6]));
mother=deval(sol,times);
states=repmat(reshape(mother,6,1,[]),1,4,1);
states(1,:,:)=states(1,:,:)+50000;
states(:,1,1)=mother(:,1);
states(:,2,513)=mother(:,513);
states(:,3,1025)=mother(:,1025);
% 碎片 4 在不同时刻满足两个条件，但从未同时满足。
states(1,4,2)=mother(1,2);
states(4,4,2)=mother(4,2)+200;
[expected,mask,first]=reference_count_cleared_debris(oe,states,times);
assert(expected==3 && isequaln(first,[0,512,1024,NaN]));
for block=[1,64,512,2048]
    [actual,actualMask,actualFirst]=count_cleared_debris(oe,states,times,block);
    assert(actual==expected && isequal(actualMask,mask) && isequaln(actualFirst,first));
end
% 单个碎片、单个时刻，避免 squeeze 引起维度歧义。
[n,~,t]=count_cleared_debris(oe,states(:,1,1),0,512);
assert(n==1 && t==0);

b=[7000e3,7004e3;.02,.02;1,1;0,0;0,0;0,2*pi];
step=[2000,1,1,1,1,pi/2];
base=struct('TopK',3,'Levels',2,'Verbose',false,'RefineThreshold',5,'ArchiveThreshold',8);
score=@(x) 8+floor((x(1)-7000e3)/1000); % 大量同分，用于检验固定排序。
expected=reference_coarse_to_fine_grid_search(score,b,step,base);
for batch=[1,7,128]
    opts=base; opts.BatchSize=batch;
    actual=coarse_to_fine_grid_search(score,b,step,opts);
    same_search(expected,actual);
    assert(actual.CacheHits>0 && actual.Evaluations<actual.GridVisits);
    assert(actual.GridVisits==expected.Evaluations);
end
opts.UseCache=false;
actual=coarse_to_fine_grid_search(score,b,step,opts);
same_search(expected,actual);
assert(actual.CacheHits==0 && actual.Evaluations==expected.Evaluations);
% 精确缓存不会把相距 1e-13 的偏心率当作同一候选。
tiny=b; tiny(1,:)=[7000e3,7000e3]; tiny(6,:)=[0,0]; tiny(2,:)=[.02,.02+1e-13];
actual=coarse_to_fine_grid_search(@(x) x(2),tiny,ones(1,6), ...
    struct('Levels',0,'Verbose',false));
assert(actual.Evaluations==2);
% 每次调用使用新缓存，不会混用评分上下文。
actual=coarse_to_fine_grid_search(@(x) -1,tiny,ones(1,6), ...
    struct('Levels',0,'Verbose',false));
assert(actual.BestScore==-1 && actual.Evaluations==2);

windowBounds=[oe',oe']; windowBounds(6,:)=[oe(6),oe(6)+.001];
windowOptions=struct('Levels',0,'Verbose',false,'ArchiveThreshold',1,'BatchSize',2);
serial=execute_debris_grid_search(states,times,windowBounds,ones(1,6), ...
    windowOptions,struct('EnableParallel',false),physics);
assert(~serial.Execution.ParallelUsed);
if ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox')
    % 清除主进程全局常量，验证工作进程使用传入 physics 初始化。
    J_2=[]; R_E=[]; mu=[];
    parallel=execute_debris_grid_search(states,times,windowBounds,ones(1,6), ...
        windowOptions,struct('EnableParallel',true,'NumWorkers',2),physics);
    assert(parallel.Execution.ParallelUsed);
    same_search(serial,parallel);
    opts=base; opts.UseParallel=true; opts.BatchSize=7; opts.MaxWorkers=2;
    actual=coarse_to_fine_grid_search(score,b,step,opts);
    same_search(expected,actual);
else
    fprintf('Parallel comparison skipped: toolbox unavailable.\n');
end
J_2=physics.J_2; R_E=physics.R_E; mu=physics.mu;
fprintf('test_grid_search_optimization passed\n');
end

function same_search(a,b)
assert(isequal(a.BestOE,b.BestOE) && a.BestScore==b.BestScore);
assert(isequal(a.Candidates,b.Candidates));
assert(isequal(a.PreferredOrbits,b.PreferredOrbits));
assert(isequal([a.History.BestScore],[b.History.BestScore]));
assert(isequal(vertcat(a.History.Step),vertcat(b.History.Step)));
assert(strcmp(a.StopReason,b.StopReason));
end
