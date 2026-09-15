function [indices, relativeTime] = debris_search_time_window(dataStart, frameCount, stepSeconds, searchStart, searchEnd)
% 根据 UTC 时刻截取碎片帧；母航天器初始根数历元为 searchStart。
% 边界必须落在现有采样点上，不进行静默取整或轨迹插值。
validateattributes(frameCount,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(stepSeconds,{'numeric'},{'scalar','positive','finite'});
dates={dataStart,searchStart,searchEnd};
for k=1:3
    assert(isdatetime(dates{k}) && isscalar(dates{k}) && ~isnat(dates{k}), ...
        '时间必须为有效的 datetime 标量。');
    assert(strcmp(dates{k}.TimeZone,'UTC'),'时间必须使用 UTC 时区。');
end
offsets=seconds([searchStart-dataStart,searchEnd-dataStart])/stepSeconds;
assert(offsets(2)>=offsets(1),'末尾时间不能早于初始时间。');
rounded=round(offsets);
assert(all(abs(offsets-rounded)<1e-7), ...
    '搜索起止时间必须与碎片数据采样点对齐（间隔 %g 秒）。',stepSeconds);
assert(rounded(1)>=0 && rounded(2)<frameCount, ...
    '指定搜索时间超出碎片轨迹覆盖范围。');
indices=(rounded(1):rounded(2))+1;
relativeTime=(0:numel(indices)-1)'*stepSeconds;
end
