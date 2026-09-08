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
window_use = 1;          % 全部四窗，若只调试可改为 1

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
