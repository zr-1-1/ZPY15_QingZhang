function finalize_mission_improvements()
% 必须在 run_mission_combination 完成后运行，避免同时写最终方案。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission');
diary(fullfile(output,'improvement_validation.log')); cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
saved=load(fullfile(output,'best_mission.mat'));
loaded=load(fullfile(output,'combination_catalog.mat'),'catalog'); catalog=loaded.catalog;
loaded=load(fullfile(output,'improved_combinations.mat'),'improvement'); improvement=loaded.improvement;
states=load(fullfile(root,'data','x_rv.mat'),'x_rv');
bestPlan=saved.bestPlan; bestVerification=saved.bestVerification;
bestScore=bestVerification.Score; bestDV=sum([bestPlan.Ships.TotalDeltaV]);
attempts=saved.attempts; baseIndex=max([attempts.Candidate]);
tight=odeset(mission_ode_options(),'RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15);
for k=1:min(8,size(improvement.Candidates,1))
    genes=improvement.Candidates(k,1:21);
    [departures,polishedScore]=polish_mission_departures(genes,catalog);
    for polished=[true,false]
        if polished, plan=assemble_mission_plan(genes,catalog,departures); plannedScore=polishedScore;
        else, plan=assemble_mission_plan(genes,catalog); plannedScore=improvement.Candidates(k,22); end
        verifiedScore=nan;
        if plan.Feasible
            verification=verify_mission_plan(plan,states.x_rv); verifiedScore=verification.Score;
            dv=sum([plan.Ships.TotalDeltaV]);
            if verification.Valid && (verifiedScore>bestScore || (verifiedScore==bestScore && dv<bestDV))
                stable=verify_mission_plan(plan,states.x_rv,tight);
                if stable.Valid && stable.Score==verifiedScore && ...
                        isequal(sort(stable.Events.DebrisID),sort(verification.Events.DebrisID))
                    bestPlan=plan; bestVerification=verification; bestScore=verifiedScore; bestDV=dv;
                else
                    plan.Failure='收紧积分容差后清除集合不一致，未作为最终候选。';
                end
            end
        end
        attempts(end+1)=struct('Candidate',baseIndex+k,'Polished',polished,'PlannedScore',plannedScore, ...
            'Feasible',plan.Feasible,'VerifiedScore',verifiedScore,'Message',plan.Failure); %#ok<AGROW>
        fprintf('补充候选 %d polished=%d planned=%d feasible=%d verified=%g best=%d。\n', ...
            k,polished,plannedScore,plan.Feasible,verifiedScore,bestScore);
    end
end
verification=verify_mission_plan(bestPlan,states.x_rv,tight);
assert(verification.Valid && verification.Score==bestScore && ...
    isequal(sort(verification.Events.DebrisID),sort(bestVerification.Events.DebrisID)),'补充方案严容差复核不一致。');
comparison=struct('OriginalScore',bestVerification.Score,'TightScore',verification.Score, ...
    'MaxAltitudeDifference_m',max(abs(verification.Altitudes-bestVerification.Altitudes),[],'all'), ...
    'SameEventTimes',isequal(verification.Events.TimeSeconds,bestVerification.Events.TimeSeconds));
saved.bestPlan=bestPlan; saved.bestVerification=verification; saved.attempts=attempts;
saved.validationComparison=comparison; saved.additionalSearch=improvement;
save(fullfile(output,'best_mission.mat'),'-struct','saved','-v7.3');
export_mission_report(bestPlan,verification,catalog,saved.searches,attempts,comparison,root);
fprintf('补充搜索后最终清除数：%d。\n',verification.Score);
end
