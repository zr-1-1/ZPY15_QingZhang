function verify_atk_adapter()
% Cross-check complete saved J2000 trajectories and the provisional correction.
global R_E
R_E=6378137;
load('data/atk_correction/corrected_plan.mat','plan','calibration');
load('data/atk_correction/dynamics_audit.mat','raw');
options=odeset(mission_ode_options(),'RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15);
files=dir('ATK_result/*J2000*.csv'); audit=zeros(3,4);
for ship=1:3
    file=files(startsWith({files.name},sprintf('Mother%d-',ship)));
    tab=readtable(fullfile(file.folder,file.name),'NumHeaderLines',2,'ReadVariableNames',false);
    time=seconds(datetime(tab{:,1},'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS')-datetime(2030,11,14,8,0,0));
    reference=tab{:,2:7}';
    burns=raw(cellfun(@(x)x.mother==ship && strcmp(x.kind,'CMCSManeuver'),raw));
    R=rotation(calibration,0); x=plan.Ships(ship).InitialState;
    x=[R*x(1:3);R*x(4:6)]; now=0; err=zeros(1,numel(time)); h=inf;
    for arc=1:7
        stop=86400;
        if arc<=6,stop=plan.Ships(ship).Pulses(arc,1);end
        sol=ode113(@(t,y)atk_calibrated_dynamics(t,y,calibration),[now,stop],x,options);
        mask=time>=now & time<=stop;
        y=deval(sol,time(mask)); ts=time(mask);
        for k=1:numel(ts),R=rotation(calibration,ts(k));y(1:3,k)=R'*y(1:3,k);end
        err(mask)=vecnorm(y(1:3,:)-reference(1:3,mask));
        h=min(h,mission_minimum_altitude(sol));
        x=sol.y(:,end);
        if arc<=6,R=rotation(calibration,stop);x(4:6)=x(4:6)+R*burns{arc}.burn;end
        now=stop;
    end
    audit(ship,:)=[ship,max(err),h,min(vecnorm(reference(1:3,:)))-R_E];
end
auditTable=array2table(audit,'VariableNames',{'Mother','MaxPositionDifference_m','AdapterMinimumAltitude_m','ATKSampleMinimumAltitude_m'});
disp(auditTable); writetable(auditTable,'data/atk_correction/full_trajectory_audit.csv');
debris=load('data/x_rv.mat','x_rv');
originalPlan=plan;
for ship=1:3
    burns=raw(cellfun(@(x)x.mother==ship && strcmp(x.kind,'CMCSManeuver'),raw));
    for k=1:6,originalPlan.Ships(ship).Pulses(k,2:4)=burns{k}.burn';end
end
originalVerification=verify_atk_plan(originalPlan,debris.x_rv,calibration,options);
fprintf('Original XML adapter replay: score %d, valid %d\n',originalVerification.Score,originalVerification.Valid);
disp(originalVerification.EligiblePerShip);
verification=verify_atk_plan(plan,debris.x_rv,calibration,options);
fprintf('Provisional corrected adapter replay: score %d, valid %d\n',verification.Score,verification.Valid);
disp(verification.MinimumAltitudes);disp(verification.EligiblePerShip);
save('data/atk_correction/adapter_verification.mat','verification','originalVerification','auditTable','-v7.3');
end
function R=rotation(c,t)
w=c.FramePolynomial*[1;t/86400]; S=[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
R=eye(3)+S+.5*S*S;
end
