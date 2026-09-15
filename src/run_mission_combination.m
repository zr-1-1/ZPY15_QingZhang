function run_mission_combination()
% 完整联合优化入口。先运行 run_seeded_mission_search，保留原有轨道库。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission');
diary(fullfile(output,'combination.log')); closeDiary=onCleanup(@()diary('off')); %#ok<NASGU>
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
sourceFiles=cell(1,5); sourceInfo=cell(1,5); libraries=cell(1,4);
sourceFiles{1}=fullfile(root,'data','preferred_orbit_library.mat');
legacy=load(sourceFiles{1}); sourceInfo{1}=dir(sourceFiles{1});
for window=1:4
    sourceFiles{window+1}=fullfile(output,sprintf('seeded_window%d.mat',window));
    sourceInfo{window+1}=dir(sourceFiles{window+1});
    loaded=load(sourceFiles{window+1},'result'); libraries{window}=loaded.result.Orbits;
end
libraries{1}=[legacy.preferredOrbitLibrary.Orbits;libraries{1}];
for window=1:4
    [~,keep]=unique(libraries{window}(:,1:6),'rows','stable'); libraries{window}=libraries{window}(keep,:);
end
states=load(fullfile(root,'data','x_rv.mat'),'x_rv');
catalogFile=fullfile(output,'combination_catalog.mat');
if isfile(catalogFile)
    saved=load(catalogFile);
    assert(isequal(saved.sourceInfo,sourceInfo),'候选库已变化，请为新一轮搜索移走旧 combination_catalog.mat。');
    catalog=saved.catalog;
else
    if isempty(gcp('nocreate'))
        cluster=parcluster('Processes'); cluster.NumWorkers=max(cluster.NumWorkers,6); cluster.NumThreads=1;
        parpool(cluster,6);
    end
    catalog=prepare_mission_catalog(libraries,states.x_rv,300);
    save(catalogFile,'catalog','sourceInfo','sourceFiles','-v7.3');
end
seeds=[1,7,23,47,89,131]; searches=cell(1,numel(seeds)); candidates=zeros(0,22);
for k=1:numel(seeds)
    searches{k}=optimize_mission_combination(catalog,seeds(k));
    candidates=[candidates;searches{k}.Genes,searches{k}.Score; ...
        searches{k}.Population,searches{k}.EliteScores]; %#ok<AGROW>
end
[~,uniqueRows]=unique(candidates(:,1:21),'rows','stable'); candidates=candidates(uniqueRows,:);
candidates=sortrows(candidates,-22);
save(fullfile(output,'combination_searches.mat'),'searches','candidates');
bestScore=-inf; bestDeltaV=inf; bestPlan=[]; bestVerification=[];
attempts=struct('Candidate',{},'Polished',{},'PlannedScore',{},'Feasible',{},'VerifiedScore',{},'Message',{});
% 对排名前 12 个不同方案复核；每个都尝试事件边界细化和原始时隙。
for k=1:min(12,size(candidates,1))
    genes=candidates(k,1:21);
    [departures,polishedScore]=polish_mission_departures(genes,catalog);
    for polished=[true,false]
        if polished, plan=assemble_mission_plan(genes,catalog,departures); plannedScore=polishedScore;
        else, plan=assemble_mission_plan(genes,catalog); plannedScore=candidates(k,22); end
        verifiedScore=nan; message=plan.Failure;
        if plan.Feasible
            verification=verify_mission_plan(plan,states.x_rv);
            verifiedScore=verification.Score;
            dv=sum([plan.Ships.TotalDeltaV]);
            if verification.Valid && (verifiedScore>bestScore || (verifiedScore==bestScore && dv<bestDeltaV))
                bestPlan=plan; bestVerification=verification; bestScore=verifiedScore; bestDeltaV=dv;
                save(fullfile(output,'best_mission_checkpoint.mat'),'bestPlan','bestVerification');
            end
        end
        attempts(end+1)=struct('Candidate',k,'Polished',polished,'PlannedScore',plannedScore, ...
            'Feasible',plan.Feasible,'VerifiedScore',verifiedScore,'Message',message); %#ok<AGROW>
        fprintf('候选 %d，时间细化=%d，计划=%d，可行=%d，复核=%g，当前最佳=%g。\n', ...
            k,polished,plannedScore,plan.Feasible,verifiedScore,bestScore);
    end
end
assert(~isempty(bestPlan),'当前候选均未通过拼接，需要扩大转移时刻或轨道组合搜索。');
% 更严容差独立复核；最终表格使用这一轮结果。
tight=odeset(mission_ode_options(),'RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15);
finalVerification=verify_mission_plan(bestPlan,states.x_rv,tight);
assert(finalVerification.Valid && finalVerification.Score==bestScore && ...
    isequal(sort(finalVerification.Events.DebrisID),sort(bestVerification.Events.DebrisID)), ...
    '更严容差验证不一致，需要检查边界事件。');
validationComparison=struct('OriginalScore',bestVerification.Score,'TightScore',finalVerification.Score, ...
    'MaxAltitudeDifference_m',max(abs(finalVerification.Altitudes-bestVerification.Altitudes),[],'all'), ...
    'SameEventTimes',isequal(finalVerification.Events.TimeSeconds,bestVerification.Events.TimeSeconds));
bestVerification=finalVerification;
save(fullfile(output,'best_mission.mat'),'bestPlan','bestVerification','searches','attempts', ...
    'validationComparison','sourceFiles','-v7.3');
export_mission_report(bestPlan,bestVerification,catalog,searches,attempts,validationComparison,root);
fprintf('最终最佳可行采样清除数：%d。\n',bestVerification.Score);
end
