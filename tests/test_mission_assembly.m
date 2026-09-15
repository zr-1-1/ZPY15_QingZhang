function test_mission_assembly()
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
initial=orb_elements2rv([7200e3,.02,98*pi/180,1.32,.7,2])';
sol=ode45(@model_rv,[0,86400],initial,mission_ode_options());
debris=repmat(reshape(deval(sol,0:86400),6,1,[]),1,8,1);
starts=[0,6,12,16]*3600; libraries=cell(1,4);
for k=1:4
    state=deval(sol,starts(k)); libraries{k}=[rv2coe(state(1:3),state(4:6)),8];
end
catalog=prepare_mission_catalog(libraries,debris,300);
plan=assemble_mission_plan(ones(1,21),catalog);
assert(plan.Feasible && plan.PlannedCoastScore==8);
assert(all([plan.Ships.TotalDeltaV]<.01));
verification=verify_mission_plan(plan,debris);
assert(verification.Valid && verification.Score==8);
assert(isequal(verification.OwnedPerShip,[8,0,0]));
fprintf('test_mission_assembly passed.\n');
end
