function test_j2_transfer()
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
initial=orb_elements2rv([7200e3,.02,98*pi/180,1.32,.7,2])';
reference=ode45(@model_rv,[0,2700],initial,mission_ode_options());
transfer=solve_j2_transfer(initial,reference.y(:,end),2700);
assert(transfer.Feasible && transfer.PositionResidual<.01);
assert(transfer.TotalDeltaV<1e-3,'同轨道自然传播应无需机动。');
target=orb_elements2rv([7300e3,.025,98.1*pi/180,1.33,.8,5])';
transfer=solve_j2_transfer(initial,target,2700);
assert(transfer.Feasible && transfer.MinAltitude>=200e3);
replay=ode45(@model_rv,[0,2700],initial+[zeros(3,1);transfer.DeltaV(:,1)],mission_ode_options());
arrival=replay.y(:,end)+[zeros(3,1);transfer.DeltaV(:,2)];
assert(norm(arrival(1:3)-target(1:3))<.01);
assert(norm(arrival(4:6)-target(4:6))<1e-7);
fprintf('test_j2_transfer passed.\n');
end
