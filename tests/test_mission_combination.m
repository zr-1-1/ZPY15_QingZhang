function test_mission_combination()
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
% 同一碎片的两次接近必须都保留，不能只用 firstTime 判断截断收益。
initial=orb_elements2rv([7200e3,.02,98*pi/180,1.32,.7,2])';
solution=ode45(@model_rv,[0,10],initial,mission_ode_options());
debris=reshape(deval(solution,0:10),6,1,11);
debris(1,1,4:8)=debris(1,1,4:8)+100000;
intervals=mission_hit_intervals(solution,debris,100,3);
assert(isequal(intervals,[1,100,102;1,108,110]));
% 合成覆盖问题：三个航天器在每窗从四个互斥单目标中选三，最优为 12。
catalog.DebrisCount=16; catalog.Orbits=repmat({1:4},1,4);
catalog.Departures={1,1,1}; catalog.FullLibraryUnion=true(16,1);
for window=1:4
    catalog.MaskSizes{window}=[4,1,1]; masks=false(16,4);
    masks((window-1)*4+(1:4),:)=eye(4)>0;
    catalog.Masks{window}=masks;
end
genes=ones(1,21); [score,~,ships]=score_mission_genes(genes,catalog);
assert(score==4 && all(sum(ships)==4),'跨航天器重复目标不得累加。');
result=optimize_mission_combination(catalog,17,struct('Population',30,'Generations',20,'Elite',3,'Immigrants',3));
assert(result.Score==12,'应达到可独立证明的合成覆盖上界。');
% 独立完整回放：两艘同轨道航天器对两个同位置碎片，全球只计两个。
solution=ode45(@model_rv,[0,86400],initial,mission_ode_options());
debris=repmat(reshape(deval(solution,0:86400),6,1,[]),1,2,1);
ship=struct('InitialState',initial,'Pulses',zeros(0,4));
plan=struct('Feasible',true,'Ships',[ship,ship], ...
    'EpochUTC',datetime(2030,11,14,8,0,0,'TimeZone','UTC'));
verified=verify_mission_plan(plan,debris);
assert(verified.Valid && verified.Score==2);
assert(isequal(verified.OwnedPerShip,[2,0]) && all(verified.Events.TimeSeconds==0));
% Propagation and debris frames may differ: rotate around the common J2 axis.
R=[0,-1,0;1,0,0;0,0,1]; Q=blkdiag(R,R);
rotated=plan;
for k=1:2,rotated.Ships(k).InitialState=Q*initial;end
frameVerified=verify_mission_plan(rotated,debris,mission_ode_options(), ...
    @mission_propagate,@(~,y)Q'*y);
assert(frameVerified.Score==2 && all(frameVerified.Events.TimeSeconds==0));
% A monotone-radius arc may have an explicitly empty event-state array.
assert(mission_minimum_altitude(struct('y',initial,'ye',[]))==norm(initial(1:3))-R_E);
% t=0 脉冲造成零长度前弧，必须同时检查脉冲后状态，且多目标维度正确。
delta=200*initial(4:6)/norm(initial(4:6));
post=initial+[zeros(3,1);delta];
solution=ode45(@model_rv,[0,86400],post,mission_ode_options());
debris=repmat(reshape(deval(solution,0:86400),6,1,[]),1,2,1);
plan.Ships=struct('InitialState',initial,'Pulses',[0,delta']);
verified=verify_mission_plan(plan,debris);
assert(verified.Score==2 && all(verified.Events.TimeSeconds==0));
fprintf('test_mission_combination passed.\n');
end
