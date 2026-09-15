function crosscheck_atk_correction()
% Validate the actual exported XML using an independent numerical integrator.
global R_E
R_E=6378137;
load('data/atk_correction/corrected_plan.mat','plan','calibration');
saved=load('data/atk_correction/adapter_verification.mat','verification','originalVerification');
raw=jsondecode(fileread('data/atk_correction/corrected_segments.json'));
if isstruct(raw),raw=num2cell(raw);end
stateError=0; impulseError=0;
for ship=1:3
    initial=raw{find(cellfun(@(x)x.mother==ship && strcmp(x.kind,'CMCSInitialState'),raw),1)};
    stateError=max(stateError,max(abs(initial.InitialState.rv-plan.Ships(ship).InitialState)));
    plan.Ships(ship).InitialState=initial.InitialState.rv;
    burns=raw(cellfun(@(x)x.mother==ship && strcmp(x.kind,'CMCSManeuver'),raw));
    for k=1:6
        impulseError=max(impulseError,max(abs(burns{k}.burn-plan.Ships(ship).Pulses(k,2:4)')));
        plan.Ships(ship).Pulses(k,2:4)=burns{k}.burn';
    end
end
assert(stateError<1e-9 && impulseError<1e-9);
reports=jsondecode(fileread('data/atk_correction/original_report_counts.json'));
for ship=1:3
    actual=sort(reports.Reports(ship).EligibleDebris);
    predicted=find(isfinite(saved.originalVerification.FirstTimes(:,ship)));
    assert(isequal(actual(:),predicted(:)),'Adapter does not reproduce original ATK clearance set');
end
debris=load('data/x_rv.mat','x_rv');
options=odeset(mission_ode_options(),'RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15);
verification=verify_atk_plan(plan,debris.x_rv,calibration,options,@ode89);
assert(verification.Valid && isequal(sort(verification.Events.DebrisID),sort(saved.verification.Events.DebrisID)));
crosscheck=struct('Score',verification.Score,'MinimumAltitudes',verification.MinimumAltitudes, ...
    'MaxStateExportError',stateError,'MaxImpulseExportError',impulseError, ...
    'MaxAltitudeDifference_m',max(abs(verification.Altitudes-saved.verification.Altitudes),[],'all'), ...
    'NativeATKVerified',false,'OriginalATKClearanceSetsReproduced',true);
disp(crosscheck);
% Decimal rounding reflects manual input/export precision, not a probabilistic guarantee.
rounding=[];
for decimals=[9,6]
    trial=plan;
    for ship=1:3,trial.Ships(ship).Pulses(:,2:4)=round(trial.Ships(ship).Pulses(:,2:4),decimals);end
    v=verify_atk_plan(trial,debris.x_rv,calibration,options,@ode89);
    rounding=[rounding;decimals,v.Score,v.Valid,v.MinimumAltitudes]; %#ok<AGROW>
end
roundingTable=array2table(rounding,'VariableNames',{'PulseDecimalPlaces','Score','Valid','Mother1Min_m','Mother2Min_m','Mother3Min_m'});
disp(roundingTable);
save('data/atk_correction/export_crosscheck.mat','crosscheck','roundingTable','-v7.3');
writetable(roundingTable,'data/atk_correction/rounding_sensitivity.csv');
writetable(verification.Events,'data/atk_correction/cleared_debris.csv');
end
