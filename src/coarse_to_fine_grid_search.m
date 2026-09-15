function result = coarse_to_fine_grid_search(scoreFcn, bounds, coarseStep, options)
% COARSE_TO_FINE_GRID_SEARCH 六维轨道网格搜索，评分越大越好。
% scoreFcn 接收 1x6 [a,e,i,Omega,omega,M]（m、rad），返回有限实数。
% bounds 为 6x2 连续展开区间；后三维按 2*pi 周期处理。
% options: TopK(3), Levels(2), RefineFactor(2), MaxEvaluations(200000),
%          Verbose(true), RefineThreshold(-Inf)。只有分数严格大于阈值才细搜。
%          ArchiveThreshold(Inf)：所有评分 >= 该值的轨道单独入库，不受 TopK 限制。
%          BatchSize(128), UseCache(true), UseParallel(false)。
% 并行模式要求调用者已建立进程池，评分函数自行初始化所需物理参数。
% result.BestOE / BestScore / Candidates / History / Evaluations。
% 不在预算耗尽后返回貌似完整的结果：超预算会在该层开始前报错。
if nargin < 4, options = struct(); end
validateattributes(bounds, {'numeric'}, {'real','finite','size',[6,2]});
validateattributes(coarseStep, {'numeric'}, {'real','finite','positive','numel',6});
assert(isa(scoreFcn,'function_handle'), 'scoreFcn 必须是函数句柄。');
assert(all(bounds(:,2) >= bounds(:,1)), '搜索下界不能大于上界。');
assert(bounds(1,1)>0 && bounds(2,1)>=0 && bounds(2,2)<1 && ...
    bounds(3,1)>=0 && bounds(3,2)<=pi, '要求 a>0、0<=e<1、0<=i<=pi。');
assert(all(bounds(4:6,2)-bounds(4:6,1)<=2*pi+1e-12), ...
    '周期角的搜索宽度不能超过 2*pi。');
defaults = struct('TopK',3,'Levels',2,'RefineFactor',2, ...
    'MaxEvaluations',200000,'Verbose',true,'RefineThreshold',-Inf, ...
    'ArchiveThreshold',Inf,'BatchSize',128,'UseCache',true,'UseParallel',false, ...
    'MaxWorkers',4,'ProgressIntervalSeconds',10);
names = fieldnames(options);
for k = 1:numel(names)
    assert(isfield(defaults,names{k}), '未知选项：%s', names{k});
    defaults.(names{k}) = options.(names{k});
end
options = defaults;
validateattributes(options.TopK, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.Levels, {'numeric'}, {'scalar','integer','nonnegative','finite'});
validateattributes(options.RefineFactor, {'numeric'}, {'scalar','integer','>=',2,'finite'});
validateattributes(options.MaxEvaluations, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.Verbose, {'logical','numeric'}, {'scalar'});
validateattributes(options.RefineThreshold, {'numeric'}, {'scalar','real','nonnan'});
validateattributes(options.ArchiveThreshold, {'numeric'}, {'scalar','real','nonnan'});
validateattributes(options.BatchSize, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.UseCache, {'logical'}, {'scalar'});
validateattributes(options.UseParallel, {'logical'}, {'scalar'});
validateattributes(options.MaxWorkers, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.ProgressIntervalSeconds, {'numeric'}, {'scalar','positive','finite'});
if options.UseParallel
    assert(~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox'), ...
        '并行评分需要 Parallel Computing Toolbox。');
    pool=gcp('nocreate');
    assert(~isempty(pool) && isa(pool,'parallel.ProcessPool'), ...
        '请先创建进程池，或设置 UseParallel=false。');
end
step = coarseStep(:);
fullCircle = false(6,1);
fullCircle(4:6) = abs(diff(bounds(4:6,:),1,2)-2*pi)<1e-12;
top = zeros(0,7); % 内部保留展开角，最后一列为评分。
evaluations = 0;
gridVisits = 0;
cache = containers.Map('KeyType','char','ValueType','double');
history = struct('Level',{},'Step',{},'Evaluations',{},'BestScore',{}, ...
    'GridVisits',{},'CacheHits',{},'WallSeconds',{});
stopReason = 'Completed';
preferred = zeros(0,7);
for level = 0:options.Levels
    if level == 0
        regions = reshape(bounds,6,2,1);
    else
        seeds = top(top(:,7)>options.RefineThreshold,:);
        if isempty(seeds)
            stopReason = 'NoEligibleSeeds';
            if options.Verbose
                fprintf('没有评分严格大于 %g 的候选，停止细搜。\n',options.RefineThreshold);
            end
            break;
        end
        previousStep = step;
        step = step/options.RefineFactor;
        regions = zeros(6,2,size(seeds,1));
        for seed = 1:size(seeds,1)
            center = seeds(seed,1:6)';
            lo = max(bounds(:,1),center-previousStep);
            hi = min(bounds(:,2),center+previousStep);
            % 整圈角允许越过原区间端点，在评分前归一化。
            lo(fullCircle) = center(fullCircle)-min(previousStep(fullCircle),pi);
            hi(fullCircle) = center(fullCircle)+min(previousStep(fullCircle),pi);
            regions(:,:,seed) = [lo,hi];
        end
    end
    grids = cell(size(regions,3),6);
    count = zeros(size(regions,3),6);
    for region = 1:size(regions,3)
        for dim = 1:6
            lo = regions(dim,1,region); hi = regions(dim,2,region);
            % 先检查数量，避免极小步长造成单轴分配失控。
            assert((hi-lo)/step(dim) <= options.MaxEvaluations, ...
                '单轴网格已超预算，请增大步长或 MaxEvaluations。');
            axis = lo:step(dim):hi;
            tolerance = 32*eps(max(1,max(abs([lo,hi]))));
            if hi-axis(end)>tolerance, axis(end+1)=hi; end
            if dim>=4 && abs(hi-lo-2*pi)<1e-12 && numel(axis)>1
                axis(end)=[]; % 整圈不重复首尾点。
            end
            grids{region,dim}=axis;
            count(region,dim)=numel(axis);
        end
    end
    planned = sum(prod(count,2));
    assert(planned <= options.MaxEvaluations-gridVisits, ...
        '第 %d 层最多需 %g 次评分，累计超预算 %g；请增大步长、减少 TopK/Levels 或提高预算。', ...
        level,planned,options.MaxEvaluations);
    startCount = evaluations;
    startVisits = gridVisits;
    levelTimer = tic;
    progressTimer = tic;
    for region = 1:size(regions,3)
        for batchStart = 0:options.BatchSize:prod(count(region,:))-1
            indices=(batchStart:min(batchStart+options.BatchSize-1,prod(count(region,:))-1))';
            remainder = indices;
            candidates = zeros(numel(indices),6);
            for dim = 1:6
                digit = mod(remainder,count(region,dim))+1;
                remainder = floor(remainder/count(region,dim));
                candidates(:,dim)=reshape(grids{region,dim}(digit),[],1);
            end
            normalized = candidates;
            normalized(:,4:6)=mod(normalized(:,4:6),2*pi);
            if options.UseCache
                % IEEE double 的精确位模式：不将近邻但不同的轨道合并。
                keys=cell(size(normalized,1),1);
                for row=1:size(normalized,1)
                    point=normalized(row,:);
                    point(point==0)=0; % 统一正负零。
                    keys{row}=reshape(num2hex(point(:))',1,[]);
                end
                [uniqueKeys,firstRows,mapping]=unique(keys,'stable');
                known=isKey(cache,uniqueKeys);
                missingRows=firstRows(~known);
                newScores=evaluate_grid_batch(scoreFcn,normalized(missingRows,:), ...
                    options.UseParallel,options.MaxWorkers);
                missingKeys=uniqueKeys(~known);
                for row=1:numel(missingKeys)
                    cache(missingKeys{row})=newScores(row);
                end
                uniqueScores=cell2mat(values(cache,uniqueKeys));
                scores=uniqueScores(mapping);
                evaluations=evaluations+numel(missingRows);
            else
                scores=evaluate_grid_batch(scoreFcn,normalized,options.UseParallel,options.MaxWorkers);
                evaluations=evaluations+size(normalized,1);
            end
            gridVisits=gridVisits+size(normalized,1);
            % 固定按原网格遍历顺序归并，不受并行任务完成先后影响。
            for row=1:size(candidates,1)
                candidate=candidates(row,:);
                score=scores(row);
                if score>=options.ArchiveThreshold
                    preferred(end+1,:)=[normalized(row,:),double(score)]; %#ok<AGROW>
                end
                % 只保存 TopK；已有优胜点在后续层继续参与比较。
                delta = abs(top(:,1:6)-candidate);
                delta(:,4:6)=abs(mod(delta(:,4:6)+pi,2*pi)-pi);
                duplicate = any(all(delta <= 1e-10*max(1,abs(candidate)),2));
                if ~duplicate && (size(top,1)<options.TopK || score>top(end,7))
                    top = sortrows([top;candidate,double(score)],-7);
                    top = top(1:min(options.TopK,size(top,1)),:);
                end
            end
            if options.Verbose && toc(progressTimer)>=options.ProgressIntervalSeconds
                fprintf('第 %d 层进度 %d/%d；实际评分 %d 次，最高分 %.12g\n', ...
                    level,gridVisits-startVisits,planned,evaluations-startCount,top(1,7));
                progressTimer=tic;
            end
        end
    end
    history(end+1)=struct('Level',level,'Step',step', ...
        'Evaluations',evaluations-startCount,'BestScore',top(1,7), ...
        'GridVisits',gridVisits-startVisits, ...
        'CacheHits',gridVisits-startVisits-evaluations+startCount, ...
        'WallSeconds',toc(levelTimer)); %#ok<AGROW>
    % 合并不同局部网格和不同层重复命中的轨道；先归一化周期角。
    if ~isempty(preferred)
        keys=round(preferred(:,1:6)./[1e-6,1e-12,1e-12,1e-12,1e-12,1e-12]);
        for dim=4:6
            keys(keys(:,dim)==round(2*pi/1e-12),dim)=0;
        end
        [~,keep]=unique(keys,'rows','stable');
        preferred=preferred(keep,:);
    end
    if options.Verbose
        fprintf('第 %d 层：评分 %d 次，缓存命中 %d 次，耗时 %.2f 秒，最高分 %.12g\n', ...
            level,evaluations-startCount,history(end).CacheHits,toc(levelTimer),top(1,7));
    end
end
top(:,4:6)=mod(top(:,4:6),2*pi);
result = struct('BestOE',top(1,1:6),'BestScore',top(1,7), ...
    'Candidates',top,'History',history,'Evaluations',evaluations, ...
    'Bounds',bounds,'Options',options,'StopReason',stopReason, ...
    'PreferredOrbits',sortrows(preferred,-7),'GridVisits',gridVisits, ...
    'CacheHits',gridVisits-evaluations);
end
