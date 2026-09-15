function intervals = mission_hit_intervals(solution,debrisStates,absoluteStart,chunkSize)
% 每一行 [碎片列编号, 首个命中秒, 最后命中秒]，保留所有连续采样命中区间。
% solution 以窗口起点为零；不能仅保留首次接近，因为截断后可再次接近。
if nargin<4, chunkSize=512; end
n=size(debrisStates,2); nt=size(debrisStates,3);
hit=false(n,nt);
for first=1:chunkSize:nt
    indices=first:min(nt,first+chunkSize-1);
    states=deval(solution,indices-1);
    relative=debrisStates(:,:,indices)-reshape(states,6,1,[]);
    hit(:,indices)=reshape(sum(relative(1:3,:,:).^2,1)<30000^2 & ...
        sum(relative(4:6,:,:).^2,1)<150^2,n,[]);
end
intervals=zeros(0,3);
for id=1:n
    changes=diff([false,hit(id,:),false]);
    starts=find(changes==1)-1+absoluteStart;
    stops=find(changes==-1)-2+absoluteStart;
    intervals=[intervals;repmat(id,numel(starts),1),starts',stops']; %#ok<AGROW>
end
end
