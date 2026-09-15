function result = reference_coarse_to_fine_grid_search(scoreFcn, bounds, coarseStep, options)
% 优化前的逐点搜索基线，用于验证同分排序、种子与轨道库一致。
% COARSE_TO_FINE_GRID_SEARCH 六维轨道网格搜索，评分越大越好。
% scoreFcn 接收 1x6 [a,e,i,Omega,omega,M]（m、rad），返回有限实数。
% bounds 为 6x2 连续展开区间；后三维按 2*pi 周期处理。
% options: TopK(3), Levels(2), RefineFactor(2), MaxEvaluations(200000),
%          Verbose(true), RefineThreshold(-Inf)。只有分数严格大于阈值才细搜。
%          ArchiveThreshold(Inf)：所有评分 >= 该值的轨道单独入库，不受 TopK 限制。
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
    'ArchiveThreshold',Inf);
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
step = coarseStep(:);
fullCircle = false(6,1);
fullCircle(4:6) = abs(diff(bounds(4:6,:),1,2)-2*pi)<1e-12;
top = zeros(0,7); % 内部保留展开角，最后一列为评分。
evaluations = 0;
history = struct('Level',{},'Step',{},'Evaluations',{},'BestScore',{});
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
    assert(planned <= options.MaxEvaluations-evaluations, ...
        '第 %d 层最多需 %g 次评分，累计超预算 %g；请增大步长、减少 TopK/Levels 或提高预算。', ...
        level,planned,options.MaxEvaluations);
    startCount = evaluations;
    for region = 1:size(regions,3)
        for index = 0:prod(count(region,:))-1
            remainder = index;
            candidate = zeros(1,6);
            for dim = 1:6
                digit = mod(remainder,count(region,dim))+1;
                remainder = floor(remainder/count(region,dim));
                candidate(dim)=grids{region,dim}(digit);
            end
            normalized = candidate;
            normalized(4:6)=mod(normalized(4:6),2*pi);
            score = scoreFcn(normalized);
            validateattributes(score, {'numeric'}, {'scalar','real','finite'}, ...
                mfilename,'scoreFcn 返回值');
            evaluations = evaluations+1;
            if score>=options.ArchiveThreshold
                preferred(end+1,:)=[normalized,double(score)]; %#ok<AGROW>
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
    end
    history(end+1)=struct('Level',level,'Step',step', ...
        'Evaluations',evaluations-startCount,'BestScore',top(1,7)); %#ok<AGROW>
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
        fprintf('第 %d 层：评分 %d 次，当前最高分 %.12g\n', ...
            level,evaluations-startCount,top(1,7));
    end
end
top(:,4:6)=mod(top(:,4:6),2*pi);
result = struct('BestOE',top(1,1:6),'BestScore',top(1,7), ...
    'Candidates',top,'History',history,'Evaluations',evaluations, ...
    'Bounds',bounds,'Options',options,'StopReason',stopReason, ...
    'PreferredOrbits',sortrows(preferred,-7));
end
