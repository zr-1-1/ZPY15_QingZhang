%% filter_high_capture.m
%  从精品数据库中筛选出【捕获数 > 7】的轨道基元，按时段分别保存
clc; clear;

%% ===== 载入精品数据库 =====
S = load('refined_all.mat');   % 含 refined_all, windows
refined_all = S.refined_all;   % 1×4 cell
windows     = S.windows;       % 4×2

ncap_min = 7;                  % ★ 筛选阈值：捕获数 > 7 即 ≥8

%% ===== 分时段筛选 =====
high_all = cell(4, 1);         % 每段筛选结果
summary  = zeros(4, 4);        % [总数, 最高, ≥8, ≥10]

fprintf('========== 捕获数 > %d 的轨道基元（分时段）==========\n\n', ncap_min);

for w = 1:4
    seg = refined_all{w};
    if isempty(seg)
        fprintf('W%d: 数据库为空\n\n', w);
        high_all{w} = [];
        continue;
    end

    ncap = [seg.n_captured];
    mask = ncap > ncap_min;                 % 严格大于 7
    high = seg(mask);

    % 按捕获数降序排序
    [~, order] = sort([high.n_captured], 'descend');
    high = high(order);

    high_all{w} = high;

    if isempty(high)
        fprintf('W%d (%.0fh-%.0fh): 无满足条件的轨道\n\n', ...
                w, windows(w,1), windows(w,2));
        continue;
    end

    ncap_h = [high.n_captured];
    summary(w, :) = [numel(high), max(ncap_h), ...
                     sum(ncap_h >= 8), sum(ncap_h >= 10)];

    fprintf('W%d (%.0fh-%.0fh): %d 条，最高 %d 颗\n', ...
            w, windows(w,1), windows(w,2), numel(high), max(ncap_h));
end

%% ===== 汇总 =====
fprintf('\n========== 汇总 ==========\n');
fprintf('%-4s %-12s %-10s %-10s %-10s %-10s\n', ...
        '时段', '窗口(h)', '总条数', '最高捕获', '≥8 条数', '≥10 条数');
for w = 1:4
    if isempty(high_all{w})
        fprintf('W%d   %-12s %-10s\n', w, ...
                sprintf('%.0f-%.0f', windows(w,1), windows(w,2)), '空');
    else
        fprintf('W%d   %-12s %-10d %-10d %-10d %-10d\n', w, ...
                sprintf('%.0f-%.0f', windows(w,1), windows(w,2)), ...
                summary(w,1), summary(w,2), summary(w,3), summary(w,4));
    end
end

%% ===== 打印每条筛选出的轨道（含捕获碎片编号）=====
fprintf('\n========== 详细列表 ==========\n');
for w = 1:4
    high = high_all{w};
    if isempty(high), continue; end
    fprintf('\n--- W%d (%.0fh-%.0fh) ---\n', w, windows(w,1), windows(w,2));
    fprintf('%-4s %-10s %-9s %-9s %-10s %-9s %-9s %-6s\n', ...
            '序号','a(km)','e','i(rad)','Ω(rad)','ω(rad)','M(rad)','捕获');
    for k = 1:numel(high)
        s = high(k);
        fprintf('#%-3d %-10.1f %-9.4f %-9.5f %-10.5f %-9.4f %-9.4f %-6d\n', ...
            k, s.oe(1)/1e3, s.oe(2), s.oe(3), s.oe(4), s.oe(5), s.oe(6), ...
            s.n_captured);
        fprintf('      碎片: %s\n', mat2str(s.captured_ids));
    end
end

%% ===== 保存 =====
save('high_capture_segments.mat', ...
     'high_all', 'windows', 'ncap_min', 'summary', '-v7.3');

fprintf('\n已保存到 high_capture_segments.mat\n');
fprintf('变量说明:\n');
fprintf('  high_all{w} : W%d 时段的筛选结果（结构数组）\n', 1:4);
fprintf('    .oe            : [a, e, i, Ω, ω, M]\n');
fprintf('    .captured_ids  : 捕获碎片编号列表\n');
fprintf('    .capture_times : 每个碎片的首次捕获时刻（秒，绝对时间）\n');
fprintf('    .n_captured    : 捕获数\n');
fprintf('    .window        : [h_start, h_end]\n');