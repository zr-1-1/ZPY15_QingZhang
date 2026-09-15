function verification = verify_mission_plan(plan,debrisStates,odeOptions,propagator,stateTransform)
% 独立回放：只读取初始位置速度与脉冲，不使用基元轨迹或拼接中间状态。
% 全部 1 s 采样点均复核，转移段也计数；每个脉冲前后均判定一次。
global R_E
if nargin<3, odeOptions=mission_ode_options(); end
if nargin<4, propagator=@mission_propagate; end
assert(plan.Feasible,'不能验证未完成的拼接。');
assert(numel(plan.Ships)<=3 && size(debrisStates,3)==86401);
n=size(debrisStates,2); ns=numel(plan.Ships);
first=inf(n,ns); distances=nan(n,ns); speeds=nan(n,ns);
altitudes=zeros(86401,ns); minimum=zeros(1,ns); replay=cell(1,ns);
for ship=1:ns
    pulses=plan.Ships(ship).Pulses;
    assert(size(pulses,1)<=6 && size(pulses,2)==4 && all(isfinite(pulses),'all'));
    assert(all(pulses(:,1)>=0 & pulses(:,1)<=86400) && all(diff(pulses(:,1))>0));
    assert(all(pulses(:,1)==round(pulses(:,1))),'当前回放脉冲时刻要求对齐 1 s 网格。');
    state=plan.Ships(ship).InitialState(:); now=0; minimum(ship)=inf;
    solutions=cell(1,size(pulses,1)+1);
    for arc=1:size(pulses,1)+1
        stop=86400;
        if arc<=size(pulses,1), stop=pulses(arc,1); end
        if stop==now
            times=now; states=state; solutions{arc}=[];
            minimum(ship)=min(minimum(ship),norm(state(1:3))-R_E);
        else
            solution=propagator([now,stop],state,odeOptions);
            solutions{arc}=solution;
            minimum(ship)=min(minimum(ship),mission_minimum_altitude(solution));
        end
        for start=now:512:stop
            times=start:min(start+511,stop);
            if stop>now, states=deval(solution,times); end
            altitudes(times+1,ship)=vecnorm(states(1:3,:))'-R_E;
            comparisonStates=states;
            if nargin>=5, comparisonStates=stateTransform(times,states); end
            relative=debrisStates(:,:,times+1)-reshape(comparisonStates,6,1,[]);
            d=reshape(sqrt(sum(relative(1:3,:,:).^2,1)),n,[]);
            v=reshape(sqrt(sum(relative(4:6,:,:).^2,1)),n,[]);
            hit=d<30000 & v<150;
            [exists,offset]=max(hit,[],2);
            candidateTimes=reshape(times(offset),[],1);
            ids=find(exists & candidateTimes<first(:,ship));
            if ~isempty(ids)
                index=sub2ind(size(d),ids,offset(ids));
                first(ids,ship)=times(offset(ids));
                distances(ids,ship)=d(index); speeds(ids,ship)=v(index);
            end
        end
        if stop>now, state=solution.y(:,end); end
        if arc<=size(pulses,1), state(4:6)=state(4:6)+pulses(arc,2:4)'; end
        now=stop;
    end
    replay{ship}=solutions;
end
[firstGlobal,owner]=min(first,[],2); ids=find(isfinite(firstGlobal));
indices=sub2ind(size(first),ids,owner(ids));
events=table(ids,owner(ids),firstGlobal(ids),plan.EpochUTC+seconds(firstGlobal(ids)), ...
    distances(indices),speeds(indices),'VariableNames', ...
    {'DebrisID','Mother','TimeSeconds','UTC','Distance_m','RelativeSpeed_mps'});
events=sortrows(events,{'TimeSeconds','Mother','DebrisID'});
verification=struct('Valid',all(minimum>=200e3),'Score',numel(ids), ...
    'Events',events,'FirstTimes',first,'MinimumAltitudes',minimum, ...
    'Altitudes',altitudes,'Replay', {replay},'SampleStepSeconds',1, ...
    'EligiblePerShip',sum(isfinite(first),1),'OwnedPerShip',accumarray(owner(ids),1,[ns,1])', ...
    'Method','Independent model_rv replay; 1 s samples and radial extrema; both sides of impulses');
end
