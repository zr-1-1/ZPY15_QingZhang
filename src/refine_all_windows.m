%% refine_all_windows.m
%  从 4 个时段的粗筛数据库中提取种子，加密 M，得到精品数据库
clc; clear;

global J_2 R_E_m mu
J_2 = 1.082627e-3;  R_E_m = 6378137;  mu = 3.986004418e14;

S = load('result_24h.mat');
S2 = load('all_segments.mat');
Dstate_t1_ra = S.Dstate_t1_ra;
Dstate_t1_va = S.Dstate_t1_va;
windows = S2.windows;            % 4×2
all_segments = S2.all_segments;  % 1×4 cell

nM = 192;                        % M 加密点数
ncap_min = 4;                    % 种子阈值

refined_all = cell(4, 1);

for w = 1:4
    seg = all_segments{w};
    if isempty(seg)
        fprintf('W%d 无种子，跳过\n', w);
        continue;
    end

    ncap = [seg.n_captured];
    seeds = seg(ncap >= ncap_min);
    fprintf('\n=== W%d 种子 %d 条 (≥%d 颗) ===\n', w, numel(seeds), ncap_min);

    if isempty(seeds)
        refined_all{w} = [];
        continue;
    end

    % 提取 (a,e,i,Ω,ω)
    base5 = zeros(numel(seeds), 5);
    for k = 1:numel(seeds)
        base5(k, :) = seeds(k).oe(1:5);
    end

    % 去重（a 保留到 1 km，其他保留到 1e-5 rad）
    key = round(base5 ./ [1e3, 1e-5, 1e-5, 1e-5, 1e-5]);
    [~, ia] = unique(key, 'rows');
    base_geom = base5(ia, :);
    fprintf('W%d 去重后独立几何轨道: %d\n', w, size(base_geom, 1));

    % 加密精筛
    tag = sprintf('R%d', w);
    refined_all{w} = refine_window( ...
        windows(w,1), windows(w,2), base_geom, ...
        Dstate_t1_ra, Dstate_t1_va, ...
        mu, J_2, R_E_m, nM, tag);
end

%% 汇总
fprintf('\n========== 精品数据库汇总 ==========\n');
fprintf('%-4s %-10s %-10s %-10s %-10s\n', ...
        '时段', '总条数', '最高捕获', '≥8 条数', '≥10 条数');
for w = 1:4
    r = refined_all{w};
    if isempty(r)
        fprintf('W%d   %-10s\n', w, '空');
        continue;
    end
    ncap = [r.n_captured];
    fprintf('W%d   %-10d %-10d %-10d %-10d\n', ...
            w, numel(r), max(ncap), sum(ncap >= 8), sum(ncap >= 10));
end

save('refined_all.mat', 'refined_all', 'windows', 'nM', 'ncap_min', '-v7.3');
fprintf('\n精品数据库已保存到 refined_all.mat\n');