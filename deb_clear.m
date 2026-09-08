%% ===================== 主程序：轨道基元数据库生成（论文 2.1 节） =====================
clear; clc; close all;

% ------------------- 0. 物理常量与任务参数 -------------------
MU  = 3.986004418e14;      % 地球引力常数 [m^3/s^2]
RE  = 6378137;             % 地球赤道半径 [m]
J2  = 1.08262668e-3;       % 带谐系数 J2（与 EGM96 一致）
clear_thresh = 8;          % 入库阈值：单条轨道在窗内清除数 >= 8
d_thresh     = 30000;      % 位置容差 30 km [m]
v_thresh     = 150;        % 速度容差 150 m/s
t_search_step = 60;        % 采样步长 60 s

% 四个时间窗口（相对于 2030-11-14 08:00:00 的秒数）
time_windows = [0,     8*3600; ...   % 08:00-16:00
                6*3600,14*3600; ...   % 14:00-22:00
                12*3600,20*3600; ...  % 20:00-04:00(+1)
                16*3600,24*3600];     % 00:00(+1)-08:00(+1)
window_use = 2;          % 全部四窗，若只调试可改为 1

% ------------------- 1. 定位并读取 ZPY15.xml -------------------
xmlCand = { fullfile(pwd,'ZPY15.xml'), ...
             };
SCEN = '';
for k = 1:numel(xmlCand)
    if exist(xmlCand{k},'file')
        SCEN = xmlCand{k};
        break;
    end
end
assert(~isempty(SCEN), '找不到 ZPY15.xml，请修改 xmlCand 列表中的路径');
fprintf('使用场景文件: %s\n', SCEN);

% 读取所有碎片（Debris*）的初始位置速度
[names, rv0] = read_debris_xml(SCEN);
nd = size(rv0,1);
fprintf('碎片总数: %d\n', nd);

% 将初始 rv 转为轨道根数 (a[m], e, i, RAAN, w, M [rad])
oe_deb = zeros(nd,6);
for d = 1:nd
    oe_deb(d,:) = rv2oe_si(rv0(d,1:3), rv0(d,4:6), MU);
end

% 打印碎片根数分布概况
alt = (oe_deb(:,1)-RE)/1000;
fprintf('碎片高度范围 %.0f ~ %.0f km\n', min(alt), max(alt));
fprintf('倾角范围 %.3f~%.3f°\n', min(oe_deb(:,3))*180/pi, max(oe_deb(:,3))*180/pi);
fprintf('RAAN 范围 %.2f~%.2f°\n', min(oe_deb(:,4))*180/pi, max(oe_deb(:,4))*180/pi);

% ------------------- 2. 逐窗口建库 -------------------
orbit_base = cell(numel(window_use),1);   % 存储每窗的基元结构体

for ww = 1:numel(window_use)
    wid = window_use(ww);
    t0w = time_windows(wid,1);
    t1w = time_windows(wid,2);
    t0_hh = floor(t0w/3600); t0_mm = floor(mod(t0w,3600)/60);
    t1_hh = floor(t1w/3600); t1_mm = floor(mod(t1w,3600)/60);
    fprintf('\n===== 时间窗 %d: %02d:%02d ~ %02d:%02d =====\n', wid, t0_hh, t0_mm, t1_hh, t1_mm);

    t_win = t0w : t_search_step : t1w;
    nt = numel(t_win);

    % ---- 2.1 预计算所有碎片在窗口内各时刻的位置/速度（向量化存储） ----
    % pos_deb_all: [nt, 3, nd], vel_deb_all: [nt, 3, nd]
    pos_deb_all = zeros(nt,3,nd);
    vel_deb_all = zeros(nt,3,nd);
    for d = 1:nd
        oeT = j2sec_prop_v(oe_deb(d,:), t_win, MU, RE, J2);   % 6 x nt
        [r, v] = oe2rv_v(oeT, MU);                             % 3 x nt
        pos_deb_all(:,:,d) = r';
        vel_deb_all(:,:,d) = v';
    end
    fprintf('碎片状态预计算完成 (nt=%d, nd=%d)\n', nt, nd);

    % ---- 2.2 确定六维网格搜索范围（聚焦碎片根数聚集区域） ----
    mu_oe = mean(oe_deb, 1);
    sig_oe = std(oe_deb, 0, 1);
    % 扩展边界：a ± 200 km, e ± 0.005, 角度 ±2°（M 全周 0~2pi）
    margin = [200e3, 0.005, 2*pi/180, 2*pi/180, 2*pi/180, pi]; % 对 M 特殊处理
    bounds = [mu_oe - margin; mu_oe + margin];
    bounds(1,6) = 0;   % M 下限固定为 0
    bounds(2,6) = 2*pi;
    % 限制 i 在合理范围 (0~pi)
    bounds(1,3) = max(0, bounds(1,3));
    bounds(2,3) = min(pi, bounds(2,3));

    % ---- 2.3 Level-1: 粗网格搜索（覆盖全部范围） ----
    step1 = [100e3, 0.002, 1*pi/180, 1*pi/180, 1*pi/180, 15*pi/180];
    grid1 = cell(1,6);
    for dim = 1:6
        grid1{dim} = bounds(1,dim) : step1(dim) : bounds(2,dim);
        if dim == 6 && grid1{dim}(end) < 2*pi
            grid1{dim} = [grid1{dim}, 2*pi];
        end
    end
    n1 = cellfun(@numel, grid1);
    fprintf('Level-1 网格规模: a=%d, e=%d, i=%d, RAAN=%d, w=%d, M=%d\n', n1(1), n1(2), n1(3), n1(4), n1(5), n1(6));
    fprintf('总网格点数 ~%d\n', prod(n1));

    % 生成网格点（使用 ndgrid）
    dims = cell(1,6);
    [dims{:}] = ndgrid(grid1{:});
    pts1 = cell2mat(cellfun(@(x) x(:), dims, 'UniformOutput', false)); % N x 6

    % 评估粗网格（利用 parfor 加速，若无并行工具箱可改 for）
    cnt1 = zeros(size(pts1,1), 1);
    fprintf('开始评估 Level-1 网格点...\n');
    tic;
    parfor k = 1:size(pts1,1)
        cnt1(k) = count_capture_vec(pts1(k,:), t_win, pos_deb_all, vel_deb_all, ...
                                    MU, RE, J2, d_thresh, v_thresh);
    end
    fprintf('Level-1 完成，耗时 %.1f s\n', toc);

    % 保留清除数 >= 5 的点作为候选（阈值稍低防止遗漏）
    keep1 = cnt1 >= 5;
    pts_cand = pts1(keep1, :);
    cnt_cand = cnt1(keep1);
    [~, ord] = sort(cnt_cand, 'descend');
    topN = min(50, sum(keep1));   % 最多取 50 个潜在最优点
    pts_top = pts_cand(ord(1:topN), :);
    fprintf('Level-1 候选 %d 个，取前 %d 个进入精细化\n', sum(keep1), topN);

% ---- 2.4 Level-2: 仅对 M 精细扫描（保持其他5维不变） ----
% 相位 M 最敏感，其他维度用粗网格的中心值即可
step_M_fine = 2*pi/180;          % M 精细步长 2°
M_fine = 0 : step_M_fine : 2*pi;
pts_fine = [];
cnt_fine = [];

fprintf('开始 Level-2 精细搜索 (仅对 M，Top %d 点)...\n', topN);
tic;
% 为加速可改用 parfor，但需要收集结果；这里用 for 保持简单
for k = 1:topN
    center = pts_top(k,:);        % 1x6
    for m = M_fine
        oe_candidate = [center(1:5), m];
        c = count_capture_vec(oe_candidate, t_win, pos_deb_all, vel_deb_all, ...
                              MU, RE, J2, d_thresh, v_thresh);
        if c >= 6                 % 候选阈值，可微调
            pts_fine = [pts_fine; oe_candidate];
            cnt_fine = [cnt_fine; c];
        end
    end
    if mod(k, 10) == 0
        fprintf('  精细搜索 %d/%d，累计候选 %d 个\n', k, topN, size(pts_fine,1));
    end
end
fprintf('Level-2 完成，耗时 %.1f s，共 %d 个候选点\n', toc, size(pts_fine,1));


 % ---- 2.5 合并并精确入库 ----
    all_pts = [pts_cand(cnt_cand>=6, :); pts_fine];
    all_pts = unique(round(all_pts, 12), 'rows');  % 数值去重
    fprintf('合并去重后待入库候选点 %d 个\n', size(all_pts,1));

    base = repmat(struct('oe0', [], 'clear_num', [], 'clear_idx', []), 0, 1);
    nbase = 0;
    tic;
    for k = 1:size(all_pts,1)
        [cnt, cidx] = count_capture_vec_idx(all_pts(k,:), t_win, pos_deb_all, vel_deb_all, ...
                                            MU, RE, J2, d_thresh, v_thresh);
        if cnt >= clear_thresh
            nbase = nbase + 1;
            base(nbase) = struct('oe0', all_pts(k,:), 'clear_num', cnt, ...
                                 'clear_idx', cidx(:)');
        end
        if mod(k, 500) == 0
            fprintf('  入库筛选 %d/%d，已入库 %d 个\n', k, size(all_pts,1), nbase);
        end
    end
    orbit_base{ww} = base;
    fprintf('时间窗%d完成：入库基元 %d 个 (>=%d)，最高清除 %d，总用时 %.1f s\n', ...
        wid, nbase, clear_thresh, max([base.clear_num]), toc);
end
   %% ===================== 4. 邻域微调（在最优基元周围精细搜索） =====================
fprintf('\n===== 开始邻域微调（寻找更高清除数）=====\n');

% ---- 4.1 找出当前最优基元 ----
max_clear_all = 0;
best_oe = [];
for ww = 1:numel(orbit_base)
    if ~isempty(orbit_base{ww})
        for k = 1:numel(orbit_base{ww})
            if orbit_base{ww}(k).clear_num > max_clear_all
                max_clear_all = orbit_base{ww}(k).clear_num;
                best_oe = orbit_base{ww}(k).oe0;
                best_win = ww;
                best_idx = k;
            end
        end
    end
end

if isempty(best_oe)
    fprintf('警告：没有找到任何基元，跳过微调\n');
else
    fprintf('当前最优：窗口%d，基元%d，清除 %d 颗\n', best_win, best_idx, max_clear_all);
    fprintf('种子轨道: a=%.2f km, e=%.4f, i=%.3f°, RAAN=%.2f°, w=%.1f°, M=%.1f°\n', ...
        best_oe(1)/1000, best_oe(2), best_oe(3)*180/pi, best_oe(4)*180/pi, ...
        best_oe(5)*180/pi, best_oe(6)*180/pi);

    % ---- 4.2 定义微调邻域范围 ----
    % 每个维度的半宽度
    refine_half = [50e3, 0.001, 0.5*pi/180, 0.5*pi/180, 0.5*pi/180, 5*pi/180];
    % 微调步长（比 Level-2 更细）
    refine_step = [10e3, 0.0002, 0.1*pi/180, 0.1*pi/180, 0.1*pi/180, 1*pi/180];

    % ---- 4.3 生成邻域网格 ----
    refine_grid = cell(1,6);
    for dim = 1:6
        lo = max(bounds(1,dim), best_oe(dim) - refine_half(dim));
        hi = min(bounds(2,dim), best_oe(dim) + refine_half(dim));
        if hi - lo < refine_step(dim)
            refine_grid{dim} = linspace(lo, hi, 3);
        else
            refine_grid{dim} = lo : refine_step(dim) : hi;
        end
        if dim == 6 && refine_grid{dim}(end) < 2*pi
            refine_grid{dim} = [refine_grid{dim}, 2*pi];
        end
    end

    % 计算邻域网格点数
    n_refine = cellfun(@numel, refine_grid);
    fprintf('邻域网格规模: a=%d, e=%d, i=%d, RAAN=%d, w=%d, M=%d\n', ...
        n_refine(1), n_refine(2), n_refine(3), n_refine(4), n_refine(5), n_refine(6));
    total_refine = prod(n_refine);
    fprintf('邻域总点数: %d\n', total_refine);

    % ---- 4.4 生成并评估邻域网格 ----
    [dR{1:6}] = ndgrid(refine_grid{:});
    pts_refine = cell2mat(cellfun(@(x) x(:), dR, 'UniformOutput', false));

    cnt_refine = zeros(size(pts_refine,1), 1);
    fprintf('开始评估邻域网格点...\n');
    tic;
    parfor k = 1:size(pts_refine,1)
        cnt_refine(k) = count_capture_vec(pts_refine(k,:), t_win, ...
                                          pos_deb_all, vel_deb_all, ...
                                          MU, RE, J2, d_thresh, v_thresh);
    end
    fprintf('邻域评估完成，耗时 %.1f s\n', toc);

    % ---- 4.5 收集邻域中的高值点（清数 >= 当前最优，或至少 >= 入库阈值） ----
    % 只保留比当前最优更高或等于的点（若等于则按碎片集合去重决定是否保留）
    better_mask = cnt_refine > max_clear_all;
    equal_mask = cnt_refine == max_clear_all & cnt_refine >= clear_thresh;
    keep_mask = better_mask | equal_mask;
    pts_new = pts_refine(keep_mask, :);
    cnt_new = cnt_refine(keep_mask);

    fprintf('邻域中找到 %d 个清除数 >= %d 的点\n', sum(keep_mask), max_clear_all);

    if ~isempty(pts_new)
        % ---- 4.6 将新点与已有基元合并去重 ----
        all_existing = [];
        for ww = 1:numel(orbit_base)
            if ~isempty(orbit_base{ww})
                for k = 1:numel(orbit_base{ww})
                    all_existing = [all_existing; orbit_base{ww}(k).oe0];
                end
            end
        end
        all_pts_new = [all_existing; pts_new];
        all_pts_new = unique(round(all_pts_new, 12), 'rows');

        % ---- 4.7 重新精确入库（只处理新增的点） ----
        % 为避免重复入库，只对 pts_new 中真正新的点进行处理
        new_added = 0;
        for k = 1:size(pts_new,1)
            % 检查是否已存在（按碎片集合判重）
            [cnt, cidx] = count_capture_vec_idx(pts_new(k,:), t_win, ...
                                                pos_deb_all, vel_deb_all, ...
                                                MU, RE, J2, d_thresh, v_thresh);
            if cnt >= clear_thresh
                % 检查碎片集合是否已存在
                is_dup = false;
                for ww = 1:numel(orbit_base)
                    if ~isempty(orbit_base{ww})
                        for j = 1:numel(orbit_base{ww})
                            if isequal(sort(cidx), sort(orbit_base{ww}(j).clear_idx))
                                is_dup = true;
                                break;
                            end
                        end
                    end
                    if is_dup, break; end
                end
                if ~is_dup
                    nbase = numel(orbit_base{best_win}) + 1;
                    orbit_base{best_win}(nbase) = struct('oe0', pts_new(k,:), ...
                                                         'clear_num', cnt, ...
                                                         'clear_idx', cidx(:)');
                    new_added = new_added + 1;
                end
            end
        end
        fprintf('邻域微调新增 %d 个基元\n', new_added);
    else
        fprintf('邻域中未发现新的高值点\n');
    end

    % ---- 4.8 重新统计最高清除数 ----
    max_final = 0;
    for ww = 1:numel(orbit_base)
        if ~isempty(orbit_base{ww})
            for k = 1:numel(orbit_base{ww})
                if orbit_base{ww}(k).clear_num > max_final
                    max_final = orbit_base{ww}(k).clear_num;
                end
            end
        end
    end
    fprintf('邻域微调完成，最终最高清除数：%d 颗\n', max_final);
end
% ------------------- 3. 保存结果 -------------------
save('orbit_base_database.mat', 'orbit_base', 'time_windows', 'window_use', ...
     'clear_thresh', 'd_thresh', 'v_thresh', '-v7.3');
fprintf('\n全部计算完成，结果已保存到 orbit_base_database.mat\n');
if ~isempty(orbit_base{1})
    fprintf('第1时间窗入库基元示例 (前5个):\n');
    for k = 1:min(5,numel(orbit_base{1}))
        fprintf('  clear=%d  debris[%s]\n', ...
            orbit_base{1}(k).clear_num, num2str(orbit_base{1}(k).clear_idx));
    end
end
%% ===================== 辅助函数定义 =====================

% ---------- 读取 XML ----------
function [names, rv] = read_debris_xml(xmlfile)
    doc = xmlread(xmlfile);
    sats = doc.getElementsByTagName('Satellite');
    names = {}; rv = zeros(0,6);
    for i = 0:sats.getLength-1
        sat = sats.item(i);
        name = char(sat.getAttribute('Name'));
        if isempty(name) || ~startsWith(name, 'Debris'), continue; end
        kids = sat.getChildNodes(); orb = [];
        for j = 0:kids.getLength-1
            ch = kids.item(j);
            if ch.getNodeType == 1 && strcmp(char(ch.getTagName), 'Orbit')
                orb = ch; break;
            end
        end
        if isempty(orb), error('Debris 缺少 Orbit: %s', name); end
        names{end+1,1} = name;
        rv(end+1,:) = [readChild(orb,'PositionX'), readChild(orb,'PositionY'), ...
                       readChild(orb,'PositionZ'), readChild(orb,'VelocityX'), ...
                       readChild(orb,'VelocityY'), readChild(orb,'VelocityZ')];
    end
end

function val = readChild(node, tag)
    kids = node.getChildNodes(); val = NaN;
    for j = 0:kids.getLength-1
        ch = kids.item(j);
        if ch.getNodeType == 1 && strcmp(char(ch.getTagName), tag)
            val = str2double(char(ch.getTextContent)); return;
        end
    end
end

% ---------- 轨道根数与 rv 相互转换 ----------
function oe = rv2oe_si(r, v, mu)
% r[m], v[m/s] -> [a,e,i,RAAN,argp,M]（a[m], 角[rad]）
    r = r(:); v = v(:);
    rn = norm(r); vn = norm(v);
    h = cross(r,v); hn = norm(h);
    inc = acos(max(-1,min(1,h(3)/hn)));
    K = [0;0;1]; nv = cross(K,h); nn = norm(nv);
    if nn < 1e-12, nv = [1;0;0]; nn = 1; end
    RAAN = acos(max(-1,min(1,nv(1)/nn)));
    if nv(2) < 0, RAAN = 2*pi - RAAN; end
    evec = (1/mu)*((vn^2 - mu/rn)*r - dot(r,v)*v);
    e = norm(evec);
    a = 1/(2/rn - vn^2/mu);
    if e < 1e-12
        argp = 0;
        u = atan2(r(3)/sin(inc), r(1)*cos(RAAN) + r(2)*sin(RAAN));
        M = wrap2pi(u);
    else
        argp = acos(max(-1,min(1,dot(nv,evec)/(nn*e))));
        if evec(3) < 0, argp = 2*pi - argp; end
        cnu = max(-1,min(1,dot(evec,r)/(e*rn)));
        nu = acos(cnu);
        if dot(r,v) < 0, nu = 2*pi - nu; end
        E = 2*atan2(sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2));
        M = wrap2pi(E - e*sin(E));
    end
    oe = [a, e, inc, wrap2pi(RAAN), wrap2pi(argp), M];
end

function [r, v] = oe2rv_si(oe, mu)
% 单点转换 [a,e,i,RAAN,argp,M] -> r[m], v[m/s]
    a = oe(1); e = oe(2); i = oe(3);
    Om = oe(4); w = oe(5); M = wrap2pi(oe(6));
    E = M;
    for k = 1:60
        dE = (E - e*sin(E) - M) / (1 - e*cos(E));
        E = E - dE;
        if abs(dE) < 1e-13, break; end
    end
    x = a*(cos(E) - e);
    y = a*sqrt(1 - e^2)*sin(E);
    fac = sqrt(mu*a)/(a*(1 - e*cos(E)));
    vx = -sin(E)*fac; vy = sqrt(1 - e^2)*cos(E)*fac;
    cO = cos(Om); sO = sin(Om); cw = cos(w); sw = sin(w);
    ci = cos(i); si = sin(i);
    R = [cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si;
         sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si;
         sw*si,              cw*si,             ci];
    r = R*[x;y;0];
    v = R*[vx;vy;0];
end

function [r, v] = oe2rv_v(oe, mu)
% 向量化版本：oe 为 6xnt 列式，返回 r/v 各 3xnt
    a = oe(1,:); e = oe(2,:); i = oe(3,:);
    Om = oe(4,:); w = oe(5,:); M = wrap2pi(oe(6,:));
    E = M;
    for k = 1:40
        dE = (E - e.*sin(E) - M) ./ (1 - e.*cos(E));
        E = E - dE;
        if max(abs(dE)) < 1e-12, break; end
    end
    cE = cos(E); sE = sin(E);
    x = a.*(cE - e);
    y = a.*sqrt(1 - e.^2).*sE;
    fac = sqrt(mu.*a)./(a.*(1 - e.*cE));
    vx = -sE.*fac;
    vy = sqrt(1 - e.^2).*cE.*fac;
    cO = cos(Om); sO = sin(Om); cw = cos(w); sw = sin(w);
    ci = cos(i); si = sin(i);
    R11 = cO.*cw - sO.*sw.*ci; R12 = -cO.*sw - sO.*cw.*ci;
    R21 = sO.*cw + cO.*sw.*ci; R22 = -sO.*sw + cO.*cw.*ci;
    R31 = sw.*si;               R32 = cw.*si;
    r = [R11.*x + R12.*y; R21.*x + R22.*y; R31.*x + R32.*y];
    v = [R11.*vx + R12.*vy; R21.*vx + R22.*vy; R31.*vx + R32.*vy];
end

% ---------- J2 一阶长期摄动（赛题式(4)） ----------
function oe_t = j2sec_prop(oe0, t, mu, RE, J2)
% 单点外推
    a = oe0(1); e = oe0(2); i = oe0(3);
    Om0 = oe0(4); w0 = oe0(5); M0 = oe0(6);
    n = sqrt(mu/a^3);
    p = a*(1 - e^2);
    fac = 1.5*J2*(RE/p)^2*n;
    Om_t = Om0 - fac*cos(i)*t;
    w_t  = w0  + 0.5*fac*(5*cos(i)^2 - 1)*t;
    M_t  = M0  + n*t + 0.75*J2*(RE/p)^2*n*sqrt(1 - e^2)*(3*cos(i)^2 - 1)*t;
    oe_t = [a, e, i, wrap2pi(Om_t), wrap2pi(w_t), wrap2pi(M_t)];
end

function oeT = j2sec_prop_v(oe0, tv, mu, RE, J2)
% 向量化外推，tv 为 1xnt，返回 6xnt
    a = oe0(1); e = oe0(2); i = oe0(3);
    Om0 = oe0(4); w0 = oe0(5); M0 = oe0(6);
    tv = tv(:).';
    n = sqrt(mu/a^3);
    p = a*(1 - e^2);
    fac = 1.5*J2*(RE/p)^2*n;
    OmT = Om0 - fac*cos(i)*tv;
    wT  = w0  + 0.5*fac*(5*cos(i)^2 - 1)*tv;
    MT  = M0  + n*tv + 0.75*J2*(RE/p)^2*n*sqrt(1 - e^2)*(3*cos(i)^2 - 1)*tv;
    oeT = [a + 0*tv; e + 0*tv; i + 0*tv; wrap2pi(OmT); wrap2pi(wT); wrap2pi(MT)];
end

% ---------- 计数函数（向量化，含 200km 高度约束） ----------
function cnt = count_capture_vec(oe0, t_win, pos_deb_all, vel_deb_all, mu, RE, J2, d_thr, v_thr)
% 统计母星沿 oe0 飞行时，在 t_win 内能清除的碎片数量
% 同时检查母星是否始终高于 200 km，若违反则直接返回 0
    oeT = j2sec_prop_v(oe0, t_win, mu, RE, J2);
    [pos_c, vel_c] = oe2rv_v(oeT, mu);
    pos_c = pos_c'; vel_c = vel_c';   % nt x 3
% 仅返回清除数量，无碎片循环
    oeT = j2sec_prop_v(oe0, t_win, mu, RE, J2);
    [pos_c, vel_c] = oe2rv_v(oeT, mu);
    pos_c = pos_c'; vel_c = vel_c';   % nt x 3
    dist_sq = sum((pos_c - pos_deb_all).^2, 2);  % nt x 1 x nd
    vel_sq  = sum((vel_c - vel_deb_all).^2, 2);
    hit = any(squeeze(dist_sq < d_thr^2 & vel_sq < v_thr^2), 1);
    cnt = sum(hit);
end

function [cnt, idx] = count_capture_vec_idx(oe0, t_win, pos_deb_all, vel_deb_all, mu, RE, J2, d_thr, v_thr)
% 返回清除数量及碎片编号
    oeT = j2sec_prop_v(oe0, t_win, mu, RE, J2);
    [pos_c, vel_c] = oe2rv_v(oeT, mu);
    pos_c = pos_c'; vel_c = vel_c';
    dist_sq = sum((pos_c - pos_deb_all).^2, 2);
    vel_sq  = sum((vel_c - vel_deb_all).^2, 2);
    hit = any(squeeze(dist_sq < d_thr^2 & vel_sq < v_thr^2), 1);
    idx = find(hit);
    cnt = numel(idx);
end

% ---------- 工具函数 ----------
function a = wrap2pi(a)
    a = mod(a, 2*pi);
end
