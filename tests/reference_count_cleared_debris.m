function [count, cleared, firstTime] = reference_count_cleared_debris(oe, debrisStates, times)
% 优化前的串行评分基线，仅供等价性测试和性能比较。
% 按同步采样点判定：距离 < 30000 m 且速度差 < 150 m/s。
% debrisStates: 6 x 碎片数 x 时间数，ECI 下 m、m/s；times 相对根数历元的秒数。
% 调用者应预先校验整份数据并初始化全局 J_2、R_E、mu。
times = times(:);
validateattributes(times,{'numeric'},{'real','finite','nonnegative','nonempty'});
assert(all(diff(times)>0),'时间必须严格递增。');
assert(size(debrisStates,1)==6 && size(debrisStates,3)==numel(times), ...
    '碎片状态必须为 6 x N x numel(times)。');
cleared = false(1,size(debrisStates,2));
firstTime = nan(size(cleared));
% 母航天器用 model_rv 的位置速度动力学积分；初始根数历元为 t=0。
initialState = orb_elements2rv(oe);
if times(end)>0
    odeOptions = odeset('RelTol',1e-9,'AbsTol', ...
        [1e-3;1e-3;1e-3;1e-6;1e-6;1e-6]);
    solution = ode45(@model_rv,[0,times(end)],initialState',odeOptions);
end
chunkSize = 512;
for first = 1:chunkSize:numel(times)
    indices = first:min(first+chunkSize-1,numel(times));
    if times(end)>0
        states = deval(solution,times(indices))';
    else
        states = initialState;
    end
    for offset = 1:numel(indices)
        t = indices(offset);
        active = find(~cleared);
        relative = debrisStates(:,active,t)-states(offset,:)';
        hit = sum(relative(1:3,:).^2,1)<30000^2 & ...
            sum(relative(4:6,:).^2,1)<150^2;
        cleared(active(hit))=true;
        firstTime(active(hit))=times(t);
        if all(cleared), break; end
    end
    if all(cleared), break; end
end
count = nnz(cleared);
end
