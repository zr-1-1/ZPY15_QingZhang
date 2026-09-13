%% search_segments.m
%  粗筛：J2 平均根数解析（快）
%  精筛：ode78 + Sat_J2ODE 完整 J2（精确，与 ATK 一致）
clc; clear; close all;

%% ===== 常量 =====
mu    = 3.986004418e14;
J_2   = sqrt(5) * 4.841653717360e-4;
R_E_m = 6378137;

t_step = 3;  t0 = 0;  t1 = 28800;
t_grid = t0 : t_step : t1;
Nt = numel(t_grid);

threshold_r = 30e3;    % 精筛
threshold_v = 150;
R_coarse    = 200e3;   % 粗筛（覆盖采样 + 模型差异）

%% ===== 载入碎片 =====
S = load('result.mat');
Dstate_t1_ra = S.Dstate_t1_ra;
Dstate_t1_va = S.Dstate_t1_va;
[~, Nt_deb, N_deb] = size(Dstate_t1_ra);
assert(Nt_deb == Nt);

%% ===== 粗筛时间网格 =====
idx_c   = 1 : 20 : Nt;   if idx_c(end) ~= Nt, idx_c = [idx_c, Nt]; end
t_c     = t_grid(idx_c);
r_deb_c = Dstate_t1_ra(:, idx_c, :);
v_deb_c = Dstate_t1_va(:, idx_c, :);
Nc      = numel(t_c);

%% ===== 网格（聚焦碎片主峰）=====
a_list     = 7120e3 : 4e3 : 7200e3;
e_list     = 0 : 0.004 : 0.020;
i_list     = 1.705 : 0.003 : 1.720;      % ★ 修正到 98° 附近
Omega_list = 1.28 : 0.008 : 1.336;
omega_list = 0 : pi/2 : 3*pi/2;
M_list     = 0 : pi/6 : 11*pi/6;

N_grid = numel(a_list)*numel(e_list)*numel(i_list)*...
         numel(Omega_list)*numel(omega_list)*numel(M_list);
fprintf('网格总数：%d\n\n', N_grid);

%% ===== 阶段 1：粗筛（J2 平均根数）=====
fprintf('=== 阶段1：解析粗筛（阈值 %.0f km）===\n', R_coarse/1e3);
tic;

coarse_oe   = zeros(N_grid, 6);
coarse_cand = cell(N_grid, 1);
n_cand = 0;
cnt    = 0;

for ia = 1:numel(a_list)
 for ie = 1:numel(e_list)
  for ii = 1:numel(i_list)
   for iO = 1:numel(Omega_list)
    for iw = 1:numel(omega_list)
     for iM = 1:numel(M_list)
        cnt = cnt + 1;
        oe0 = [a_list(ia), e_list(ie), i_list(ii), ...
               Omega_list(iO), omega_list(iw), M_list(iM)];

        [r_sc_c, ~] = j2_mean_propagate_vec(oe0, t_c, mu, J_2, R_E_m);

        diff  = r_sc_c - r_deb_c;
        d2    = sum(diff.^2, 1);
        dmin2 = squeeze(min(d2, [], 2));
        hit_jd = find(dmin2 < R_coarse^2);

        if ~isempty(hit_jd)
            n_cand = n_cand + 1;
            coarse_oe(n_cand, :) = oe0;
            coarse_cand{n_cand}   = hit_jd(:)';
        end
     end
    end
   end
  end
 end
end

coarse_oe   = coarse_oe(1:n_cand, :);
coarse_cand = coarse_cand(1:n_cand);

fprintf('粗筛完成：%d / %d 条有候选 (%.2f%%)，耗时 %.1f s\n\n', ...
        n_cand, N_grid, 100*n_cand/N_grid, toc);

if n_cand == 0, error('粗筛无命中，检查网格'); end

%% ===== 阶段 2：精筛（ode78 + Sat_J2ODE）=====
fprintf('=== 阶段2：数值精筛（ode78，判据 30 km / 150 m/s）===\n');
tic;

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% 用临时数组存，最后转 struct（parfor 友好）
seg_oe    = zeros(n_cand, 6);
seg_cap   = cell(n_cand, 1);
seg_time  = cell(n_cand, 1);
seg_ncap  = zeros(n_cand, 1);
n_seg     = 0;

parfor k = 1:n_cand        % 若无 Parallel Toolbox，改为 for
    oe0 = coarse_oe(k, :);
    [Sat_r, Sat_v] = orb_elements2rv(oe0, mu);
    state0 = [Sat_r; Sat_v];

    % ★ 完整 J2 数值积分，与 ATK 一致
    [~, Y] = ode78(@Sat_J2ODE, t_grid, state0, opts);
    r_sc = Y(:, 1:3)';
    v_sc = Y(:, 4:6)';

    captured = [];
    ctimes   = [];

    for jd = coarse_cand{k}
        r_db = Dstate_t1_ra(:, :, jd);
        v_db = Dstate_t1_va(:, :, jd);

        dr = sqrt(sum((r_sc - r_db).^2, 1));
        dv = sqrt(sum((v_sc - v_db).^2, 1));

        hit_t = find(dr < threshold_r & dv < threshold_v);
        if ~isempty(hit_t)
            captured(end+1) = jd;            %#ok<AGROW>
            ctimes(end+1)   = hit_t(1);      %#ok<AGROW>
        end
    end

    if ~isempty(captured)
        seg_oe(k, :)  = oe0;
        seg_cap{k}    = captured;
        seg_time{k}   = t_grid(ctimes);
        seg_ncap(k)   = numel(captured);
    end

    if mod(k, 500) == 0
        fprintf('  精筛进度：%d / %d\n', k, n_cand);
    end
end

fprintf('精筛完成，耗时 %.1f s\n\n', toc);

%% ===== 组装 segments 结构数组 =====
hit = find(seg_ncap > 0);
n_seg = numel(hit);
segments = struct('oe', {}, 'captured_ids', {}, ...
                  'capture_times', {}, 'n_captured', {});
for kk = 1:n_seg
    k = hit(kk);
    segments(kk) = struct( ...
        'oe',            seg_oe(k, :), ...
        'captured_ids',  seg_cap{k}, ...
        'capture_times', seg_time{k}, ...
        'n_captured',    seg_ncap(k));
end

fprintf('数值精筛后：%d 条轨道满足 30 km / 150 m/s\n\n', n_seg);

%% ===== 统计 + Top 榜 =====
if n_seg == 0
    fprintf('无捕获。可能原因：\n');
    fprintf('  1) i_list 单位不是 rad（检查 1.705~1.720）\n');
    fprintf('  2) 粗筛阈值 R_coarse 太小漏检\n');
    fprintf('  3) 网格分辨率不够\n');
    return;
end

ncap = [segments.n_captured];
fprintf('=== 捕获数分布 ===\n');
for n = 1:max(ncap)
    c = sum(ncap == n);
    if c > 0, fprintf('  捕获 %2d 个：%4d 条\n', n, c); end
end

[~, order] = sort(ncap, 'descend');
fprintf('\n=== Top 20 轨道段 ===\n');
fprintf('%-4s %-10s %-9s %-9s %-10s %-9s %-9s %-6s\n', ...
        '排名','a(km)','e','i(rad)','Ω(rad)','ω(rad)','M(rad)','捕获');
for k = 1:min(20, n_seg)
    idx = order(k); s = segments(idx);
    fprintf('#%-3d %-10.1f %-9.4f %-9.5f %-10.5f %-9.4f %-9.4f %-6d\n', ...
        k, s.oe(1)/1e3, s.oe(2), s.oe(3), s.oe(4), s.oe(5), s.oe(6), s.n_captured);
end

save('segments_nopulse.mat', 'segments', 'R_coarse', ...
     'threshold_r', 'threshold_v', 't_grid', '-v7.3');
fprintf('\n已保存到 segments_nopulse.mat\n');