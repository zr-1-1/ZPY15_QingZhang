function crosscheck_mission_export()
% 检查导出精度，并用另一高阶积分器 ode89 独立回放最终 CSV 方案。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission');
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
saved=load(fullfile(output,'best_mission.mat'));
initial=readtable(fullfile(output,'initial_orbits.csv'));
pulses=readtable(fullfile(output,'maneuvers.csv'));
plan=saved.bestPlan; maxStateError=0; maxPulseError=0; maxPositionError=0; maxVelocityError=0;
for ship=1:3
    state=initial{initial.Mother==ship,{'x_m','y_m','z_m','vx_mps','vy_mps','vz_mps'}}';
    impulse=pulses{pulses.Mother==ship,{'TimeSeconds','DeltaVx_mps','DeltaVy_mps','DeltaVz_mps'}};
    maxStateError=max(maxStateError,max(abs(state-plan.Ships(ship).InitialState)));
    maxPositionError=max(maxPositionError,max(abs(state(1:3)-plan.Ships(ship).InitialState(1:3))));
    maxVelocityError=max(maxVelocityError,max(abs(state(4:6)-plan.Ships(ship).InitialState(4:6))));
    maxPulseError=max(maxPulseError,max(abs(impulse-plan.Ships(ship).Pulses),[],'all'));
    plan.Ships(ship).InitialState=state; plan.Ships(ship).Pulses=impulse;
end
assert(maxStateError<1e-7 && maxPulseError<1e-8,'CSV 精度不足，需要提高导出精度。');
states=load(fullfile(root,'data','x_rv.mat'),'x_rv');
options=odeset(mission_ode_options(),'RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15);
verification=verify_mission_plan(plan,states.x_rv,options,@independent_propagator);
assert(verification.Valid && verification.Score==saved.bestVerification.Score && ...
    isequal(sort(verification.Events.DebrisID),sort(saved.bestVerification.Events.DebrisID)), ...
    'CSV + ode89 与最终 MAT + ode113 的清除集合不同。');
crosscheck=struct('Solver','ode89','Score',verification.Score, ...
    'MaxExportStateError',maxStateError,'MaxExportPulseError',maxPulseError, ...
    'MaxExportPositionError_m',maxPositionError,'MaxExportVelocityError_mps',maxVelocityError, ...
    'MaxAltitudeDifference_m',max(abs(verification.Altitudes-saved.bestVerification.Altitudes),[],'all'), ...
    'MinimumAltitudes',verification.MinimumAltitudes,'Events',verification.Events);
save(fullfile(output,'export_crosscheck.mat'),'crosscheck');
saved.validationComparison.Crosscheck=crosscheck;
save(fullfile(output,'best_mission.mat'),'-struct','saved','-v7.3');
catalog=load(fullfile(output,'combination_catalog.mat'),'catalog');
export_mission_report(saved.bestPlan,saved.bestVerification,catalog.catalog,saved.searches, ...
    saved.attempts,saved.validationComparison,root);
fprintf('CSV + ode89 独立交叉验证通过：%d 个碎片。\n',verification.Score);
end

function solution = independent_propagator(times,state,options)
solution=ode89(@model_rv,times,state,options);
end
