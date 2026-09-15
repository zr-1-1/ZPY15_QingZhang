function correct_atk_plan()
% Generate a starting correction using the empirical adapter. Native ATK follows.
load('data/atk_correction/calibration.mat','calibration');
load('data/atk_correction/dynamics_audit.mat','raw');
old=load('data/mission/best_mission.mat','bestPlan','bestVerification');
plan=old.bestPlan; targets=zeros(6,3,3);
opts=odeset(mission_ode_options(),'RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15);
for ship=1:3
    initial=raw{find(cellfun(@(x)x.mother==ship && strcmp(x.kind,'CMCSInitialState'),raw),1)};
    state=initial.InitialState.rv; plan.Ships(ship).InitialState=state;
    R=rotation(calibration,0); state=[R*state(1:3);R*state(4:6)]; now=0;
    for gap=1:3
        departure=plan.Ships(ship).Pulses(2*gap-1,1); arrival=departure+2700;
        coast=propagate([now,departure],state,calibration,opts); state=coast.y(:,end);
        target=old.bestVerification.Replay{ship}{2*gap+1}.y(:,1);
        R=rotation(calibration,arrival); target=[R*target(1:3);R*target(4:6)];
        targets(:,gap,ship)=target;
        R=rotation(calibration,departure);
        velocity=state(4:6)+R*old.bestPlan.Ships(ship).Pulses(2*gap-1,2:4)';
        for iteration=1:12
            sol=propagate([departure,arrival],[state(1:3);velocity],calibration,opts);
            e=sol.y(1:3,end)-target(1:3);
            if norm(e)<1e-5, break; end
            jac=zeros(3);
            for j=1:3
                q=velocity; q(j)=q(j)+.01;
                trial=propagate([departure,arrival],[state(1:3);q],calibration,opts);
                jac(:,j)=(trial.y(1:3,end)-sol.y(1:3,end))/.01;
            end
            velocity=velocity-jac\e;
        end
        assert(norm(e)<1e-4,'Adapter shooting failed');
        dv1=R'*(velocity-state(4:6)); R=rotation(calibration,arrival);
        dv2=R'*(target(4:6)-sol.y(4:6,end));
        plan.Ships(ship).Pulses(2*gap-1,2:4)=dv1';
        plan.Ships(ship).Pulses(2*gap,2:4)=dv2';
        state=[sol.y(1:3,end);target(4:6)]; now=arrival;
    end
end
save('data/atk_correction/corrected_plan.mat','plan','targets','calibration');
p=[];
for ship=1:3,p=[p;repmat(ship,6,1),(1:6)',plan.Ships(ship).Pulses];end %#ok<AGROW>
writetable(array2table(p,'VariableNames',{'Mother','Pulse','TimeSeconds','DeltaVx_mps','DeltaVy_mps','DeltaVz_mps'}), ...
    'data/atk_correction/maneuvers.csv');
fprintf('Adapter correction generated; native ATK verification pending.\n');
end
function R=rotation(c,t)
w=c.FramePolynomial*[1;t/86400]; R=expm([0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0]);
end
function sol=propagate(t,x,c,o)
sol=ode113(@(t,y)atk_calibrated_dynamics(t,y,c),t,x,o);
end
