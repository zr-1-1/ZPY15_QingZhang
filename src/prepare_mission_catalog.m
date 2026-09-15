function catalog = prepare_mission_catalog(orbitLibraries,debrisStates,departureStep)
% 将四个轨道库转换为带全部命中区间、传播解和截断掩码的组合优化数据。
if nargin<3, departureStep=300; end
global J_2 R_E mu
validateattributes(departureStep,{'double'},{'scalar','positive','integer'});
assert(numel(orbitLibraries)==4);
catalog.Windows=[0,8;6,14;12,20;16,24]*3600;
catalog.TransferSeconds=2700; catalog.DepartureStep=departureStep;
catalog.DebrisCount=size(debrisStates,2);
catalog.Departures=cell(1,3);
for gap=1:3
    lo=catalog.Windows(gap+1,1); hi=catalog.Windows(gap,2)-catalog.TransferSeconds;
    catalog.Departures{gap}=unique([lo:departureStep:hi,hi]);
end
catalog.Orbits=cell(1,4); catalog.Masks=cell(1,4); catalog.MaskSizes=cell(1,4);
for window=1:4
    start=catalog.Windows(window,1); stop=catalog.Windows(window,2);
    incoming=start; outgoing=stop;
    if window>1, incoming=catalog.Departures{window-1}+catalog.TransferSeconds; end
    if window<4, outgoing=catalog.Departures{window}; end
    points=orbitLibraries{window}(:,1:6); entries=cell(1,size(points,1));
    context=struct('States',debrisStates(:,:,start+1:stop+1),'Start',start,'Stop',stop, ...
        'Physics',struct('J_2',J_2,'R_E',R_E,'mu',mu));
    pool=[];
    if ~isempty(ver('parallel')), pool=gcp('nocreate'); end
    if ~isempty(pool) && isa(pool,'parallel.ProcessPool')
        data=parallel.pool.Constant(context);
        for first=1:64:size(points,1)
            last=min(first+63,size(points,1)); batch=points(first:last,:); part=cell(1,size(batch,1));
            parfor (row=1:size(batch,1),min(6,pool.NumWorkers))
                part{row}=make_entry(batch(row,:),first+row-1,data.Value);
            end
            entries(first:last)=part;
            fprintf('复核窗口 %d：%d/%d。\n',window,last,size(points,1));
        end
        delete(data);
    else
        for row=1:size(points,1)
            entries{row}=make_entry(points(row,:),row,context);
        end
    end
    orbits=[entries{:}];
    assert(~isempty(orbits),'窗口 %d 无符合 >=8 门槛的轨道，需要补充搜索。',window);
    dims=[numel(orbits),numel(incoming),numel(outgoing)];
    masks=false(catalog.DebrisCount,prod(dims));
    for orbit=1:dims(1)
        intervals=orbits(orbit).Intervals;
        for arrival=1:dims(2)
            for departure=1:dims(3)
                retained=intervals(:,3)>=incoming(arrival) & intervals(:,2)<=outgoing(departure);
                index=sub2ind(dims,orbit,arrival,departure);
                masks(unique(intervals(retained,1)),index)=true;
            end
        end
    end
    catalog.Orbits{window}=orbits; catalog.Masks{window}=masks; catalog.MaskSizes{window}=dims;
    fprintf('组合库窗口 %d：%d 条，单段最佳 %d，%d 种截断状态。\n', ...
        window,numel(orbits),max([orbits.Score]),prod(dims));
end
catalog.FullLibraryUnion=false(catalog.DebrisCount,1);
for window=1:4
    for orbit=1:numel(catalog.Orbits{window})
        catalog.FullLibraryUnion(catalog.Orbits{window}(orbit).Intervals(:,1))=true;
    end
end
end

function entry = make_entry(oe,row,context)
global J_2 R_E mu
J_2=context.Physics.J_2; R_E=context.Physics.R_E; mu=context.Physics.mu;
entry=[];
solution=ode45(@model_rv,[0,context.Stop-context.Start],orb_elements2rv(oe)',mission_ode_options());
altitude=mission_minimum_altitude(solution);
if altitude<200e3, return; end
intervals=mission_hit_intervals(solution,context.States,context.Start);
count=numel(unique(intervals(:,1)));
if count<8, return; end
entry=struct('OE',oe,'SourceRow',row,'Solution',solution,'Intervals',intervals, ...
    'Score',count,'MinAltitude',altitude);
end
