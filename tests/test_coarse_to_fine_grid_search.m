function test_coarse_to_fine_grid_search
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
epoch=datetime(2030,11,14,8,0,0,'TimeZone','UTC');
[indices,times]=debris_search_time_window(epoch,11,60, ...
    epoch+minutes(2),epoch+minutes(5));
assert(isequal(indices,3:6) && isequal(times,[0;60;120;180]));
for badWindow={[-60,60],[0,660],[61,120],[180,120]}
    failed=false;
    try
        window=badWindow{1};
        debris_search_time_window(epoch,11,60, ...
            epoch+seconds(window(1)),epoch+seconds(window(2)));
    catch
        failed=true;
    end
    assert(failed,'非法时间窗口应报错。');
end
opts = struct('TopK',1,'Levels',2,'Verbose',false);
b = [7000e3,7004e3; .02,.02; 1,1; 0,0; 0,0; 0,0];
r = coarse_to_fine_grid_search(@(x) -((x(1)-7001500)/1000)^2, ...
    b,[2000,1,1,1,1,1],opts);
assert(abs(r.BestOE(1)-7001500)<1e-6 && r.BestScore==0);
assert(all(diff([r.History.BestScore])>=0));
% 整圈去重与越过零度细搜。
b(1,:) = 7000e3; b(6,:)=[0,2*pi];
r = coarse_to_fine_grid_search(@(x) -abs(mod(x(6)+pi/8+pi,2*pi)-pi), ...
    b,[1000,1,1,1,1,pi/2],opts);
assert(r.History(1).Evaluations==4);
assert(abs(r.BestOE(6)-15*pi/8)<1e-10);
% 展开的部分角区间覆盖上界，不越过原有范围。
b(6,:)=[345,445]*pi/180;
r = coarse_to_fine_grid_search(@(x) -abs(x(6)-85*pi/180), ...
    b,[1000,1,1,1,1,60*pi/180],opts);
assert(abs(r.BestOE(6)-85*pi/180)<1e-10);
assert(r.History(1).Evaluations==3);
opts.MaxEvaluations=1;
failed=false;
try
    coarse_to_fine_grid_search(@(x) 0,b,[1000,1,1,1,1,pi/4],opts);
catch
    failed=true;
end
assert(failed,'预算不足应明确报错。');
opts.MaxEvaluations=200000;
opts.RefineThreshold=5;
r=coarse_to_fine_grid_search(@(x) 5,b,[1000,1,1,1,1,pi/4],opts);
assert(numel(r.History)==1 && strcmp(r.StopReason,'NoEligibleSeeds'));
r=coarse_to_fine_grid_search(@(x) 6,b,[1000,1,1,1,1,pi/4],opts);
assert(numel(r.History)==3);
% 入库取 >= 门槛，不受 TopK 限制；各层重复轨道去重。
archiveBounds=[7000e3,7002e3;.02,.02;1,1;0,0;0,0;0,0];
archiveOpts=struct('TopK',1,'Levels',1,'RefineThreshold',5, ...
    'ArchiveThreshold',8,'Verbose',false);
r=coarse_to_fine_grid_search(@(x) 8,archiveBounds, ...
    [1000,1,1,1,1,1],archiveOpts);
assert(size(r.PreferredOrbits,1)==4 && all(r.PreferredOrbits(:,7)==8));
assert(size(r.Candidates,1)==1);
r=coarse_to_fine_grid_search(@(x) 7,archiveBounds, ...
    [1000,1,1,1,1,1],archiveOpts);
assert(isempty(r.PreferredOrbits) && numel(r.History)==2);
% 实际评分：重复命中只计一次；距离/速度必须在同一时刻都达标。
global J_2 R_E mu
J_2=1.0826e-3; R_E=6378137; mu=3.986004418e14;
oe=[7000e3,.02,1,0,0,0];
times=[0;1]; states=zeros(6,5,2);
initial=orb_elements2rv(oe)';
sol=ode45(@model_rv,[0,1],initial,odeset('RelTol',1e-9, ...
    'AbsTol',[1e-3;1e-3;1e-3;1e-6;1e-6;1e-6]));
for t=1:2
    sat=deval(sol,times(t));
    states(:,:,t)=repmat(sat,1,5);
    states(1,2,t)=sat(1)+31000; % 距离不合格
    states(4,3,t)=sat(4)+151; % 速度不合格
    states(1,4,t)=sat(1)+31000;
    if t==1
        states(4,5,t)=sat(4)+151;
    else
        states(1,5,t)=sat(1)+31000;
    end
end
[n,mask,first]=count_cleared_debris(oe,states,times);
assert(n==1 && isequal(mask,[true,false,false,false,false]));
assert(first(1)==0 && all(isnan(first(2:end))));
edgeStates=repmat(initial,1,2);
edgeStates(1,1)=initial(1)+30000;
edgeStates(4,2)=initial(4)+150;
assert(count_cleared_debris(oe,edgeStates,0)==0,'边界必须严格小于。');
fprintf('test_coarse_to_fine_grid_search passed\n');
end
