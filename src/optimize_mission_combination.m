function result = optimize_mission_combination(catalog,seed,options)
% 带精英保留、单点交叉、随机换段/换时隙和移民的联合最大覆盖搜索。
% 这里只评价截断后的基元；真实转移可行性及转移段收益由后续全轨迹复核。
if nargin<2, seed=1; end
if nargin<3, options=struct(); end
defaults=struct('Population',100,'Generations',400,'Elite',8,'MutationRate',.15, ...
    'Immigrants',10,'PolishPasses',4);
fields=fieldnames(options);
for k=1:numel(fields)
    assert(isfield(defaults,fields{k}),'未知选项 %s',fields{k});
    defaults.(fields{k})=options.(fields{k});
end
o=defaults;
assert(o.Population>o.Elite+o.Immigrants && o.Elite>=1);
previous=rng; cleanup=onCleanup(@()rng(previous)); %#ok<NASGU>
rng(seed,'twister');
limits=[repelem(cellfun(@numel,catalog.Orbits),3),repelem(cellfun(@numel,catalog.Departures),3)];
population=floor(rand(o.Population,21).*limits)+1;
history=zeros(o.Generations,2); timer=tic;
for generation=1:o.Generations
    scores=zeros(o.Population,1);
    for row=1:o.Population, scores(row)=score_mission_genes(population(row,:),catalog); end
    [scores,order]=sort(scores,'descend'); population=population(order,:);
    history(generation,:)=[scores(1),mean(scores)];
    if mod(generation,50)==0 || generation==1
        fprintf('组合 seed=%d generation=%d best=%d mean=%.2f elapsed=%.1fs\n', ...
            seed,generation,scores(1),mean(scores),toc(timer));
    end
    if generation==o.Generations, break; end
    next=population;
    for row=o.Elite+1:o.Population-o.Immigrants
        candidates=randi(o.Population,2,3); parents=min(candidates,[],2);
        split=randi(20);
        child=[population(parents(1),1:split),population(parents(2),split+1:end)];
        mutate=rand(1,21)<o.MutationRate;
        child(mutate)=floor(rand(1,nnz(mutate)).*limits(mutate))+1;
        next(row,:)=child;
    end
    next(end-o.Immigrants+1:end,:)=floor(rand(o.Immigrants,21).*limits)+1;
    population=next;
end
% 对优胜者逐坐标穷举，保证在已尝试邻域中不放过严格改进。
best=population(1,:); bestScore=scores(1);
for pass=1:o.PolishPasses
    improved=false;
    for gene=randperm(21)
        for value=1:limits(gene)
            trial=best; trial(gene)=value;
            score=score_mission_genes(trial,catalog);
            if score>bestScore, best=trial; bestScore=score; improved=true; end
        end
    end
    if ~improved, break; end
end
result=struct('Genes',best,'Score',bestScore,'Seed',seed,'Options',o, ...
    'History',history,'WallSeconds',toc(timer),'Population',population(1:o.Elite,:), ...
    'EliteScores',scores(1:o.Elite),'LibraryUnionUpperBound',nnz(catalog.FullLibraryUnion), ...
    'Optimality','Best found; no global optimality certificate');
end
