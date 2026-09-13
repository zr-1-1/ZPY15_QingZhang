E_0 = Debris_oe(1,:)
E_t = OE_scl_ptb(E_0,1)
[r,v] = orb_elements2rv(E_t,mu)


%% diag_w2.m
clc; clear;
load('result_24h.mat');
fprintf('=== 1. 数据尺寸 ===\n');
fprintf('Dstate_t1_ra: %s\n', mat2str(size(Dstate_t1_ra)));
fprintf('应该是 3×28801×345 (24h, 步长3s)\n');
fprintf('若是 3×9601×345，说明只生成了 8h 数据\n\n');

fprintf('=== 2. 碎片在几个时刻的位置范围 ===\n');
% t=0h, 6h, 12h, 18h 对应索引 1, 7201, 14401, 21601
for h = [0, 6, 12, 18]
    idx = h*3600/3 + 1;
    if idx > size(Dstate_t1_ra, 2)
        fprintf('  t=%2dh: 索引 %d 越界！\n', h, idx);
        continue;
    end
    r_all = squeeze(Dstate_t1_ra(:, idx, :));   % 3×345
    r_norm = sqrt(sum(r_all.^2, 1));
    fprintf('  t=%2dh: |r| ∈ [%.0f, %.0f] km, 平均 %.0f km\n', ...
            h, min(r_norm)/1e3, max(r_norm)/1e3, mean(r_norm)/1e3);
end

fprintf('\n=== 3. W2 时段碎片索引检查 ===\n');
idx_w = 7201 : 16801;
fprintf('W2 索引范围: [%d, %d], 数据最大索引: %d\n', ...
        idx_w(1), idx_w(end), size(Dstate_t1_ra, 2));
if idx_w(end) > size(Dstate_t1_ra, 2)
    fprintf('!!! 越界 !!!\n');
end

fprintf('\n=== 4. 用 W1 最优轨道在 W1/W2 时段分别测最小距离 ===\n');
S1 = load('segments_nopulse.mat');
ncap = [S1.segments.n_captured];
[~, best] = max(ncap);
oe0 = S1.segments(best).oe;
fprintf('使用 oe = [a=%.1f km, e=%.4f, i=%.4f, Ω=%.4f, ω=%.4f, M=%.4f]\n', ...
    oe0(1)/1e3, oe0(2), oe0(3), oe0(4), oe0(5), oe0(6));

mu = 3.986004418e14;
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% W1 时段
t1 = 0:3:28800; idx1 = 1:9601;
r1 = Dstate_t1_ra(:, idx1, :);
[Sat_r, Sat_v] = orb_elements2rv(oe0, mu);
[~, Y] = ode78(@Sat_J2ODE, t1, [Sat_r; Sat_v], opts);
r_sc = Y(:, 1:3)';
dmin1 = inf;
for jd = 1:345
    dr = sqrt(sum((r_sc - r1(:,:,jd)).^2, 1));
    dmin1 = min(dmin1, min(dr));
end
fprintf('W1 时段：最小距离 = %.1f km\n', dmin1/1e3);

% W2 时段
t2 = 21600:3:50400; idx2 = 7201:16801;
if idx2(end) <= size(Dstate_t1_ra, 2)
    r2 = Dstate_t1_ra(:, idx2, :);
    [Sat_r, Sat_v] = orb_elements2rv(oe0, mu);
    [~, Y] = ode78(@Sat_J2ODE, t2, [Sat_r; Sat_v], opts);
    r_sc = Y(:, 1:3)';
    dmin2 = inf;
    for jd = 1:345
        dr = sqrt(sum((r_sc - r2(:,:,jd)).^2, 1));
        dmin2 = min(dmin2, min(dr));
    end
    fprintf('W2 时段：最小距离 = %.1f km\n', dmin2/1e3);
end        
       
%% verify_fix.m
clc; clear;
load('result_24h.mat');
load('segments_nopulse.mat');
mu = 3.986004418e14;

% 用 W1 最优轨道，验证它在 W2 时段用【相对时间】能否接近
ncap = [segments.n_captured];
[~, best] = max(ncap);
oe0 = segments(best).oe;
fprintf('测试轨道: a=%.1f km, e=%.4f, i=%.4f, Ω=%.4f, ω=%.4f, M=%.4f\n', ...
    oe0(1)/1e3, oe0(2), oe0(3), oe0(4), oe0(5), oe0(6));

% W2 窗口：h_start=6, h_end=14
h_s = 6;  h_e = 14;
t_rel = 0:3:(h_e-h_s)*3600;                    % 相对时间
idx_w = (h_s*3600/3 + 1) : (h_e*3600/3 + 1);   % 绝对索引
r_deb = Dstate_t1_ra(:, idx_w, :);
v_deb = Dstate_t1_va(:, idx_w, :);

% 关键：oe0 解释为“窗口起点（h=6）时刻的根数”
[Sat_r, Sat_v] = orb_elements2rv(oe0, mu);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, Y] = ode78(@Sat_J2ODE, t_rel, [Sat_r; Sat_v], opts);
r_sc = Y(:, 1:3)';  v_sc = Y(:, 4:6)';

dmin = inf; v_at_dmin = inf;
for jd = 1:345
    dr = sqrt(sum((r_sc - r_deb(:,:,jd)).^2, 1));
    dv = sqrt(sum((v_sc - v_deb(:,:,jd)).^2, 1));
    [dm, im] = min(dr);
    if dm < dmin
        dmin = dm;  v_at_dmin = dv(im);
    end
end
fprintf('W2 时段：最小距离 = %.1f km，该点相对速度 = %.1f m/s\n', ...
    dmin/1e3, v_at_dmin);
if dmin < 30e3 && v_at_dmin < 150
    fprintf('✅ 满足判据，修复有效\n');
else
    fprintf('⚠️  仍不满足，继续诊断\n');
end