function plan = assemble_mission_plan(genes,catalog,departures)
% 从真实连续状态逐段拼接，不将交会位置强制重置为目标基元位置。
orbit=reshape(genes(1:12),3,4); slots=reshape(genes(13:21),3,3);
if nargin<3
    departures=zeros(3);
    for gap=1:3, departures(:,gap)=catalog.Departures{gap}(slots(:,gap)); end
end
assert(isequal(size(departures),[3,3]) && all(isfinite(departures),'all') && ...
    all(departures==round(departures),'all'));
for gap=1:3
    assert(all(departures(:,gap)>=catalog.Windows(gap+1,1)) && ...
        all(departures(:,gap)+catalog.TransferSeconds<=catalog.Windows(gap,2)));
end
plan=struct('Feasible',true,'Failure','','Genes',genes,'Ships',[]);
plan.EpochUTC=datetime(2030,11,14,8,0,0,'TimeZone','UTC');
for ship=1:3
    initialOE=catalog.Orbits{1}(orbit(ship,1)).OE;
    initial=orb_elements2rv(initialOE)'; state=initial; now=0;
    pulses=zeros(6,4); residuals=zeros(1,3); branches=zeros(1,3);
    minimum=inf;
    for gap=1:3
        departure=departures(ship,gap);
        arrival=departure+catalog.TransferSeconds;
        coast=mission_propagate([now,departure],state,mission_ode_options());
        minimum=min(minimum,mission_minimum_altitude(coast));
        state=coast.y(:,end);
        target=deval(catalog.Orbits{gap+1}(orbit(ship,gap+1)).Solution, ...
            arrival-catalog.Windows(gap+1,1));
        transfer=solve_j2_transfer(state,target,catalog.TransferSeconds);
        if ~transfer.Feasible
            plan.Feasible=false;
            plan.Failure=sprintf('航天器 %d，第 %d 次转移不可行，出发 t=%g s。',ship,gap,departure);
            return;
        end
        pulses(2*gap-1,:)=[departure,transfer.DeltaV(:,1)'];
        pulses(2*gap,:)=[arrival,transfer.DeltaV(:,2)'];
        state=transfer.Solution.y(:,end)+[zeros(3,1);transfer.DeltaV(:,2)];
        now=arrival;
        residuals(gap)=transfer.PositionResidual; branches(gap)=transfer.Branch;
        minimum=min(minimum,transfer.MinAltitude);
    end
    coast=mission_propagate([now,86400],state,mission_ode_options());
    minimum=min(minimum,mission_minimum_altitude(coast));
    if minimum<200e3
        plan.Feasible=false; plan.Failure=sprintf('航天器 %d 滑行段高度低于 200 km。',ship); return;
    end
    entry=struct('InitialOE',initialOE,'InitialState',initial,'Pulses',pulses, ...
        'OrbitIndices',orbit(ship,:),'SlotIndices',slots(ship,:), ...
        'PositionResiduals',residuals,'Branches',branches, ...
        'MinAltitude',minimum,'TotalDeltaV',sum(vecnorm(pulses(:,2:4),2,2)));
    plan.Ships=[plan.Ships,entry];
end
plan.Departures=departures;
plan.PlannedCoastScore=score_mission_schedule(orbit,departures,catalog);
end
