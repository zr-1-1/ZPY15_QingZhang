%% refine_seeds.m
%  以现有高捕获数轨道为种子，加密 M 相位，寻找更高捕获数的轨道段
clc; clear;

%% ===== 常量 =====
mu    = 3.986004418e14;
t_step = 3;  t0 = 0;  t1 = 28800;
t_grid = t0 : t_step : t1;
Nt = numel(t_grid);
threshold_r = 30e3;
threshold_v = 150;

%% ===== 载入 =====
S1 = load('segments_nopulse.mat', 'segments');
S2 = load('result.mat', 'Dstate_t1_ra', 'Dstate_t1_va');
segments    = S1.segments;
Dstate_t1_ra = S2.Dstate_t1_ra;
Dstate_t1_va = S2.Dstate_t1_va;

%% ===== 提取种子 =====
ncap = [segments.n_captured];
ncap_min = 4;                          % 阈值：捕获 ≥ 4 的才当种子
seed_idx = find(ncap >= ncap_min);
fprintf('种子数（捕获 ≥ %d）：%d\n', ncap_min, numel(seed_idx));

% 提取 (a, e, i, Ω, ω)
base_list = zeros(numel(seed_idx), 5);
for k = 1:numel(seed_idx)
    base_list(k, :) = segments(seed_idx(k)).oe(1:5);
end

% 按 (a,e,i,Ω,ω) 去重（同一几何轨道不同 M 算重复）
[~, ia] = unique(round(base_list, 5), 'rows');
base_list = base_list(ia, :);
n_seed = size(base_list, 1);
fprintf('去重后独立几何轨道：%d 条\n', n_seed);

%% ===== M 加密网格 =====
M_list = 0 : pi/48 : 2*pi - pi/48;     % 96 个点，步长 3.75°
nM = numel(M_list);
total = n_seed * nM;
fprintf('加密网格总数：%d\n\n', total);

% 组装所有候选 oe
all_oe = zeros(total, 6);
k = 0;
for i = 1:n_seed
    for j = 1:nM
        k = k + 1;
        all_oe(k, :) = [base_list(i,:), M_list(j)];
    end
end

%% ===== 精筛（数值积分）=====
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

tmp_cap  = cell(total, 1);
tmp_time = cell(total, 1);
tmp_n    = zeros(total, 1);

fprintf('精筛开始...\n');
tic;

parfor k = 1:total            % 若无 Parallel Toolbox，改为 for
    oe0 = all_oe(k, :);
    [Sat_r, Sat_v] = orb_elements2rv(oe0, mu);
    [~, Y] = ode78(@Sat_J2ODE, t_grid, [Sat_r; Sat_v], opts);
    r_sc = Y(:, 1:3)';
    v_sc = Y(:, 4:6)';

    captured = [];  ctimes = [];
    for jd = 1:345
        r_db = Dstate_t1_ra(:, :, jd);
        v_db = Dstate_t1_va(:, :, jd);
        dr = sqrt(sum((r_sc - r_db).^2, 1));
        dv = sqrt(sum((v_sc - v_db).^2, 1));
        hit_t = find(dr < threshold_r & dv < threshold_v, 1);
        if ~isempty(hit_t)
            captured(end+1) = jd;         %#ok<AGROW>
            ctimes(end+1)   = hit_t(1);   %#ok<AGROW>
        end
    end

    tmp_cap{k}  = captured;
    tmp_time{k} = t_grid(ctimes);
    tmp_n(k)    = numel(captured);

    if mod(k, 1000) == 0
        fprintf('  %d / %d (%.1f s)\n', k, total, toc);
    end
end

fprintf('精筛完成，耗时 %.1f s\n\n', toc);

%% ===== 组装结果 =====
hit = find(tmp_n > 0);
n_res = numel(hit);
results = struct('oe', {}, 'captured_ids', {}, ...
                 'capture_times', {}, 'n_captured', {});
for kk = 1:n_res
    k = hit(kk);
    results(kk) = struct( ...
        'oe',            all_oe(k, :), ...
        'captured_ids',  tmp_cap{k}, ...
        'capture_times', tmp_time{k}, ...
        'n_captured',    tmp_n(k));
end

%% ===== 统计 + Top 榜 =====
ncap2 = [results.n_captured];
fprintf('=== 捕获数分布 ===\n');
for n = 1:max(ncap2)
    c = sum(ncap2 == n);
    if c > 0, fprintf('  捕获 %2d 个：%4d 条\n', n, c); end
end

[~, order] = sort(ncap2, 'descend');
fprintf('\n=== Top 20 ===\n');
fprintf('%-4s %-10s %-9s %-9s %-10s %-9s %-9s %-6s\n', ...
        '排名','a(km)','e','i(rad)','Ω(rad)','ω(rad)','M(rad)','捕获');
for k = 1:min(20, n_res)
    s = results(order(k));
    fprintf('#%-3d %-10.1f %-9.4f %-9.5f %-10.5f %-9.4f %-9.4f %-6d\n', ...
        k, s.oe(1)/1e3, s.oe(2), s.oe(3), s.oe(4), s.oe(5), s.oe(6), s.n_captured);
end

%% ===== 保存 =====
save('refined_segments.mat', 'results', 'base_list', '-v7.3');
fprintf('\n已保存到 refined_segments.mat\n');