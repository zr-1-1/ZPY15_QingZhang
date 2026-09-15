function run_reference_combination_baseline()
% 比较基线：只使用已复核论文种子及原有第一窗口库，不含新增局部网格。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission');
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
loaded=load(fullfile(output,'reference_seed_diagnostic.mat'),'seeds'); libraries=cell(1,4);
for window=1:4
    selected=loaded.seeds([loaded.seeds.Window]==window & [loaded.seeds.Score]>=8);
    libraries{window}=[vertcat(selected.OE),[selected.Score]'];
end
legacy=load(fullfile(root,'data','preferred_orbit_library.mat'));
libraries{1}=[libraries{1};legacy.preferredOrbitLibrary.Orbits];
states=load(fullfile(root,'data','x_rv.mat'),'x_rv');
catalog=prepare_mission_catalog(libraries,states.x_rv,300);
searches=cell(1,3); bestScore=-inf; bestPlan=[];
for k=1:3
    searches{k}=optimize_mission_combination(catalog,k);
    [departures,~]=polish_mission_departures(searches{k}.Genes,catalog);
    plan=assemble_mission_plan(searches{k}.Genes,catalog,departures);
    if plan.Feasible
        verification=verify_mission_plan(plan,states.x_rv);
        if verification.Valid && verification.Score>bestScore
            bestScore=verification.Score; bestPlan=plan; bestVerification=verification;
        end
    end
end
assert(~isempty(bestPlan),'基线候选拼接不可行。');
save(fullfile(output,'reference_combination_baseline.mat'),'searches','bestScore','bestPlan','bestVerification');
fprintf('参考种子组合基线，独立复核清除数 %d。\n',bestScore);
end
