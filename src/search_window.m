function segments = search_window(h_start, h_end, ...
                                  Dstate_t1_ra, Dstate_t1_va, ...
                                  mu, J_2, R_E_m, tag)
% 时间约定：oe0 = 窗口起点 h_start 时刻的轨道根数
% 碎片索引、粗筛输出、精筛输出，全部对齐到绝对时刻

t_step = 3;
t_win_abs = (h_start*3600) : t_step : (h_end*3600);
Nt_w = numel(t_win_abs);

% 碎片窗口（绝对索引）
idx_w = (h_start*3600/t_step + 1) : (h_end*3600/t_step + 1);
r_deb_w = Dstate_t1_ra(:, idx_w, :);
v_deb_w = Dstate_t1_va(:, idx_w, :);

% 粗筛时间网格（60s）
idx_c = 1:20:Nt_w;
if idx_c(end) ~= Nt_w, idx_c = [idx_c, Nt_w]; end
r_deb_c  = r_deb_w(:, idx_c, :);
t_c_abs  = t_win_abs(idx_c);
t_c_rel  = t_c_abs - h_start*3600;             % ★ 相对时间

threshold_r = 30e3;
threshold_v = 150;
R_coarse    = 200e3;

a_list     = 7120e3 : 4e3 : 7200e3;
e_list     = 0 : 0.004 : 0.020;
i_list     = 1.705 : 0.003 : 1.720;
Omega_list = 1.28 : 0.008 : 1.336;
omega_list = 0 : pi/2 : 3*pi/2;
M_list     = 0 : pi/6 : 11*pi/6;               % 先用 12 点验证修复

N_grid = numel(a_list)*numel(e_list)*numel(i_list)*...
         numel(Omega_list)*numel(omega_list)*numel(M_list);
fprintf('[%s] 网格 %d，窗口 %.0fh-%.0fh\n', tag, N_grid, h_start, h_end);

% ---------- 粗筛 ----------
coarse_oe   = zeros(N_grid, 6);
coarse_cand = cell(N_grid, 1);
n_cand = 0; cnt = 0;
tic;
for ia = 1:numel(a_list)
 for ie = 1:numel(e_list)
  for ii = 1:numel(i_list)
   for iO = 1:numel(Omega_list)
    for iw = 1:numel(omega_list)
     for iM = 1:numel(M_list)
        cnt = cnt + 1;
        oe0 = [a_list(ia), e_list(ie), i_list(ii), ...
               Omega_list(iO), omega_list(iw), M_list(iM)];

        % ★ 关键：传相对时间，oe0 是 h_start 时刻的根数
        [r_sc_c, ~] = j2_mean_propagate_vec(oe0, t_c_rel, mu, J_2, R_E_m);

        d2 = sum((r_sc_c - r_deb_c).^2, 1);
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
fprintf('[%s] 粗筛 %d / %d，耗时 %.1f s\n', tag, n_cand, N_grid, toc);

% ---------- 精筛 ----------
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
seg_oe = zeros(n_cand, 6); seg_cap = cell(n_cand, 1);
seg_time = cell(n_cand, 1); seg_ncap = zeros(n_cand, 1);
tic;
parfor k = 1:n_cand
    oe0 = coarse_oe(k, :);
    [Sat_r, Sat_v] = orb_elements2rv(oe0, mu);

    % ★ 关键：用绝对时间 tspan，初值在 tspan(1)=h_start*3600 时刻
    [~, Y] = ode78(@Sat_J2ODE, t_win_abs, [Sat_r; Sat_v], opts);
    r_sc = Y(:, 1:3)';  v_sc = Y(:, 4:6)';

    captured = [];  ctimes_abs = [];
    for jd = coarse_cand{k}
        dr = sqrt(sum((r_sc - r_deb_w(:,:,jd)).^2, 1));
        dv = sqrt(sum((v_sc - v_deb_w(:,:,jd)).^2, 1));
        hit_t = find(dr < threshold_r & dv < threshold_v, 1);
        if ~isempty(hit_t)
            captured(end+1) = jd;                      %#ok<AGROW>
            ctimes_abs(end+1) = t_win_abs(hit_t(1));   %#ok<AGROW>
        end
    end
    if ~isempty(captured)
        seg_oe(k,:) = oe0;  seg_cap{k} = captured;
        seg_time{k} = ctimes_abs;  seg_ncap(k) = numel(captured);
    end
end
fprintf('[%s] 精筛耗时 %.1f s\n', tag, toc);

hit = find(seg_ncap > 0);
n_seg = numel(hit);
segments = struct('oe', {}, 'captured_ids', {}, ...
                  'capture_times', {}, 'n_captured', {}, 'window', {});
for kk = 1:n_seg
    k = hit(kk);
    segments(kk) = struct( ...
        'oe',            seg_oe(k, :), ...
        'captured_ids',  seg_cap{k}, ...
        'capture_times', seg_time{k}, ...
        'n_captured',    seg_ncap(k), ...
        'window',        [h_start, h_end]);
end

ncap = [segments.n_captured];
fprintf('[%s] 完成：%d 条，最高捕获 %d，≥4 的 %d 条\n\n', ...
    tag, n_seg, max([0, ncap]), sum(ncap >= 4));
end