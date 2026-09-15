function run_seeded_mission_search()
% 论文种子补充：两种末列角解释分别实测，严格 >5 才作局部细化，>=8 入库。
% 分块局部网格依次搜索 a/e/M 和 i/Omega/omega，保留不同覆盖集合。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission');
diary(fullfile(output,'seeded_search.log')); closeDiary=onCleanup(@()diary('off')); %#ok<NASGU>
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
physics=struct('J_2',J_2,'R_E',R_E,'mu',mu);
diagnostic=load(fullfile(output,'reference_seed_diagnostic.mat'),'seeds');
loaded=load(fullfile(root,'data','x_rv.mat'),'x_rv');
legacy=load(fullfile(root,'data','preferred_orbit_library.mat'));
windows=[0,8;6,14;12,20;16,24]*3600;
pool=gcp('nocreate');
if isempty(pool)
    cluster=parcluster('Processes'); cluster.NumWorkers=max(cluster.NumWorkers,6); cluster.NumThreads=1;
    pool=parpool(cluster,6);
end
assert(isa(pool,'parallel.ProcessPool'),'需要进程池。');
for window=1:4
    filename=fullfile(output,sprintf('seeded_window%d.mat',window));
    if isfile(filename), fprintf('保留已完成种子库 %d。\n',window); continue; end
    timer=tic;
    data=parallel.pool.Constant(struct('States',loaded.x_rv(:,:,windows(window,1)+1:windows(window,2)+1), ...
        'Times',(0:diff(windows(window,:)))','Physics',physics));
    selected=diagnostic.seeds([diagnostic.seeds.Window]==window);
    points=vertcat(selected.OE);
    if window==1, points=[points;legacy.preferredOrbitLibrary.Orbits(:,1:6)]; end
    [scores,masks]=evaluate(points,data);
    allPoints=points; allScores=scores; allMasks=masks;
    history=zeros(0,4);
    [x,y,z]=ndgrid(-1:1); offsets=[x(:),y(:),z(:)];
    for round=1:4
        dimensions=[1,2,6]; steps=[2000,.001,.25*pi/180];
        if mod(round,2)==0, dimensions=[3,4,5]; steps=[.01,.02,1]*pi/180; end
        steps=steps/2^floor((round-1)/2);
        seeds=diverse_seeds(allPoints,allScores,allMasks,12);
        assert(~isempty(seeds),'窗口 %d 没有 >5 的种子，不能进入细搜。',window);
        trials=zeros(size(seeds,1)*size(offsets,1),6);
        for k=1:size(seeds,1)
            rows=(k-1)*size(offsets,1)+(1:size(offsets,1));
            trials(rows,:)=repmat(seeds(k,:),numel(rows),1);
            trials(rows,dimensions)=trials(rows,dimensions)+offsets.*steps;
        end
        trials(:,4:6)=mod(trials(:,4:6),2*pi);
        trials=unique(trials,'rows','stable');
        % 缓存只合并数值完全相同的点；不舍入合并近邻。
        trials=trials(~ismember(trials,allPoints,'rows'),:);
        trials=trials(trials(:,1)>0 & trials(:,2)>=0 & trials(:,2)<1 & trials(:,3)>=0 & trials(:,3)<=pi,:);
        [newScores,newMasks]=evaluate(trials,data);
        allPoints=[allPoints;trials]; allScores=[allScores;newScores]; allMasks=[allMasks;newMasks]; %#ok<AGROW>
        history(end+1,:)=[round,size(trials,1),max(allScores),toc(timer)]; %#ok<AGROW>
        fprintf('种子窗口 %d round=%d evaluations=%d best=%d archive=%d elapsed=%.1fs\n', ...
            window,round,numel(allScores),max(allScores),nnz(allScores>=8),toc(timer));
    end
    result=struct('Orbits',[allPoints(allScores>=8,:),allScores(allScores>=8)], ...
        'Masks',allMasks(allScores>=8,:),'AllPoints',allPoints,'AllScores',allScores, ...
        'History',history,'RefineThreshold',5,'ArchiveThreshold',8, ...
        'WindowSeconds',windows(window,:),'Physics',physics,'WallSeconds',toc(timer), ...
        'Method','Reference seeds + two-block local grids; four rounds; 12 coverage-diverse seeds');
    save(filename,'result'); delete(data);
end
end

function seeds = diverse_seeds(points,scores,masks,limit)
eligible=find(scores>5); [~,order]=sort(scores(eligible),'descend'); eligible=eligible(order);
if isempty(eligible), seeds=zeros(0,6); return; end
% 每种完整覆盖集合先保留一个，余量再由高分点补齐。
[~,uniqueCoverage]=unique(masks(eligible,:),'rows','stable');
chosen=eligible(uniqueCoverage); chosen=chosen(1:min(limit,numel(chosen)));
remaining=setdiff(eligible,chosen,'stable');
chosen=[chosen;remaining(1:min(limit-numel(chosen),numel(remaining)))];
seeds=points(chosen,:);
end

function [scores,masks] = evaluate(points,data)
scores=zeros(size(points,1),1); masks=false(size(points,1),size(data.Value.States,2));
for first=1:128:size(points,1)
    last=min(first+127,size(points,1));
    batch=points(first:last,:); counts=zeros(size(batch,1),1); hits=false(size(batch,1),size(masks,2));
    parfor (row=1:size(batch,1),6)
        [counts(row),hits(row,:)]=evaluate_one(batch(row,:),data.Value);
    end
    scores(first:last)=counts; masks(first:last,:)=hits;
    fprintf('局部评分 %d/%d，当前批最佳 %d。\n',last,size(points,1),max(counts));
end
end

function [score,mask] = evaluate_one(oe,context)
global J_2 R_E mu
J_2=context.Physics.J_2; R_E=context.Physics.R_E; mu=context.Physics.mu;
[score,mask]=count_cleared_debris(oe,context.States,context.Times,512);
end
