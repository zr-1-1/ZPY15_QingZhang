function refined = refine_window(h_start, h_end, base_geom, ...
                                 Dstate_t1_ra, Dstate_t1_va, ...
                                 mu, J_2, R_E_m, nM, tag)
% REFINE_WINDOW 给定种子几何轨道(a,e,i,Ω,ω)，加密 M 相位，精筛
%   base_geom: N×5, 每行 [a, e, i, Ω, ω]
%   nM      : M 采样点数（推荐 192 或 288）
%   refined : 输出结构数组（只保留捕获数 > 0 的）

t_step = 3;
t_win_abs = (h_start*3600) : t_step : (h_end*3600);
Nt_w = numel(t_win_abs);

idx_w = (h_start*3600/t_step + 1) : (h_end*3600/t_step + 1);
r_deb_w = Dstate_t1_ra(:, idx_w, :);
v_deb_w = Dstate_t1_va(:, idx_w, :);

M_list = 0 : 2*pi/nM : 2*pi - 2*pi/nM;

N_geom = size(base_geom, 1);
total = N_geom * nM;
fprintf('[%s] 几何轨道 %d，M 加密 %d 点，总网格 %d\n', ...
        tag, N_geom, nM, total);

% 组装所有 oe
all_oe = zeros(total, 6);
k = 0;
for i = 1:N_geom
    for j = 1:nM
        k = k + 1;
        all_oe(k, :) = [base_geom(i, :), M_list(j)];
    end
end

threshold_r = 30e3;
threshold_v = 150;
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

tmp_cap  = cell(total, 1);
tmp_time = cell(total, 1);
tmp_n    = zeros(total, 1);

fprintf('[%s] 精筛开始...\n', tag);
tic;
parfor k = 1:total
    oe0 = all_oe(k, :);
    [Sat_r, Sat_v] = orb_elements2rv(oe0, mu);
    [~, Y] = ode78(@Sat_J2ODE, t_win_abs, [Sat_r; Sat_v], opts);
    r_sc = Y(:, 1:3)';  v_sc = Y(:, 4:6)';

    captured = [];  ctimes_abs = [];
    for jd = 1:345
        dr = sqrt(sum((r_sc - r_deb_w(:,:,jd)).^2, 1));
        dv = sqrt(sum((v_sc - v_deb_w(:,:,jd)).^2, 1));
        hit_t = find(dr < threshold_r & dv < threshold_v, 1);
        if ~isempty(hit_t)
            captured(end+1) = jd;                    %#ok<AGROW>
            ctimes_abs(end+1) = t_win_abs(hit_t(1)); %#ok<AGROW>
        end
    end
    tmp_cap{k}  = captured;
    tmp_time{k} = ctimes_abs;
    tmp_n(k)    = numel(captured);
end
fprintf('[%s] 精筛耗时 %.1f s\n', tag, toc);

hit = find(tmp_n > 0);
n_res = numel(hit);
refined = struct('oe', {}, 'captured_ids', {}, ...
                 'capture_times', {}, 'n_captured', {}, 'window', {});
for kk = 1:n_res
    k = hit(kk);
    refined(kk) = struct( ...
        'oe',            all_oe(k, :), ...
        'captured_ids',  tmp_cap{k}, ...
        'capture_times', tmp_time{k}, ...
        'n_captured',    tmp_n(k), ...
        'window',        [h_start, h_end]);
end

ncap = [refined.n_captured];
fprintf('[%s] 完成：%d 条，最高捕获 %d，≥8 的 %d 条\n\n', ...
        tag, n_res, max([0, ncap]), sum(ncap >= 8));
end