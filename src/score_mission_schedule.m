function [score,mask] = score_mission_schedule(orbitIndices,departures,catalog)
% 任意整数秒转移时刻的精确采样区间覆盖评分，用于事件边界时间细化。
mask=false(catalog.DebrisCount,1);
for ship=1:3
    incoming=[0,departures(ship,:)+catalog.TransferSeconds];
    outgoing=[departures(ship,:),86400];
    for window=1:4
        intervals=catalog.Orbits{window}(orbitIndices(ship,window)).Intervals;
        retained=intervals(:,3)>=incoming(window) & intervals(:,2)<=outgoing(window);
        mask(intervals(retained,1))=true;
    end
end
score=nnz(mask);
end
