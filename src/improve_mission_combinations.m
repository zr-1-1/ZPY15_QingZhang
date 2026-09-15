function improve_mission_combinations()
% 在已有组合库上进行低成本中性变异搜索；相同清除数也接受，跨越平台。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission');
loaded=load(fullfile(output,'combination_catalog.mat'),'catalog'); catalog=loaded.catalog;
loaded=load(fullfile(output,'combination_searches.mat'),'candidates');
candidates=loaded.candidates; best=candidates(1,1:21); bestScore=candidates(1,22);
limits=[repelem(cellfun(@numel,catalog.Orbits),3),repelem(cellfun(@numel,catalog.Departures),3)];
previous=rng; cleanup=onCleanup(@()rng(previous)); %#ok<NASGU>
rng(20260913,'twister'); restarts=24; steps=200; timer=tic; history=zeros(restarts,3);
for restart=1:restarts
    current=best;
    if mod(restart,3)==0
        current=candidates(randi(size(candidates,1)),1:21);
    end
    changed=randperm(21,randi([2,6]));
    current(changed)=floor(rand(1,numel(changed)).*limits(changed))+1;
    for step=1:steps
        gene=randi(21); scores=zeros(1,limits(gene));
        trial=current;
        for value=1:limits(gene)
            trial(gene)=value; scores(value)=score_mission_genes(trial,catalog);
        end
        top=max(scores); tied=find(scores==top);
        current(gene)=tied(randi(numel(tied)));
        if top>=bestScore
            best=current; bestScore=top;
            if mod(step,10)==0 || isempty(candidates) || top>max(candidates(:,22))
                candidates=[candidates;canonical(current),top]; %#ok<AGROW>
            end
        end
    end
    history(restart,:)=[restart,bestScore,toc(timer)];
    fprintf('中性变异 restart=%d best=%d elapsed=%.1fs\n',restart,bestScore,toc(timer));
end
candidates=[candidates;canonical(best),bestScore];
for k=1:size(candidates,1), candidates(k,1:21)=canonical(candidates(k,1:21)); end
[~,keep]=unique(candidates(:,1:21),'rows','stable'); candidates=sortrows(candidates(keep,:),-22);
improvement=struct('Candidates',candidates,'History',history,'Seed',20260913, ...
    'Restarts',restarts,'StepsPerRestart',steps,'WallSeconds',toc(timer),'BestCoastScore',bestScore);
save(fullfile(output,'improved_combinations.mat'),'improvement');
end

function genes = canonical(genes)
rows=sortrows([reshape(genes(1:12),3,4),reshape(genes(13:21),3,3)]);
genes=[reshape(rows(:,1:4),1,12),reshape(rows(:,5:7),1,9)];
end
