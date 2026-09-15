function [departures,score] = polish_mission_departures(genes,catalog)
% 固定轨道组合，穷举每次转移涉及的命中区间边界，保留严格改进。
% 整数秒评分只会在这些边界变化，不需要按秒传播或盲目密化全局网格。
orbits=reshape(genes(1:12),3,4); slots=reshape(genes(13:21),3,3);
departures=zeros(3);
for gap=1:3, departures(:,gap)=catalog.Departures{gap}(slots(:,gap)); end
score=score_mission_schedule(orbits,departures,catalog);
for pass=1:4
    improved=false;
    for ship=1:3
        for gap=1:3
            left=catalog.Orbits{gap}(orbits(ship,gap)).Intervals;
            right=catalog.Orbits{gap+1}(orbits(ship,gap+1)).Intervals;
            lo=catalog.Windows(gap+1,1); hi=catalog.Windows(gap,2)-catalog.TransferSeconds;
            candidates=unique([lo;hi;departures(ship,gap);left(:,2); ...
                right(:,3)-catalog.TransferSeconds;right(:,3)-catalog.TransferSeconds+1]);
            candidates=candidates(candidates>=lo & candidates<=hi);
            for time=candidates'
                trial=departures; trial(ship,gap)=time;
                value=score_mission_schedule(orbits,trial,catalog);
                if value>score, departures=trial; score=value; improved=true; end
            end
        end
    end
    if ~improved, break; end
end
end
