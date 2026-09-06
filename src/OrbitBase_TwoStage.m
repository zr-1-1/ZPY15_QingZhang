%% ==============================================================
clear; clc; close all;

%% ===================== 0. 常量（全 SI） =====================
MU  = 3.986004418e14;     % m^3/s^2
RE  = 6378137;            % m
C20 = -4.841653717360e-04;% EGM96 完全归一化 C20
J2  = -sqrt(2*2+1)*C20;   % 与队长 Data_Extract.m 一致(约 1.0826e-3)

clear_thresh = 8;   % 入库阈值：单条轨道在窗内清除数 >= 8
d_thresh     = 30000;  % 30 km
v_thresh     = 150;    % 150 m/s

% 时间窗（相对 2030-11-14 08:00 的秒数）
time_windows = [0,   8*3600; ...
                6*3600,  14*3600; ...
                12*3600, 20*3600; ...
                16*3600, 24*3600];
window_use = 3;         % 快速验证先只跑第 1 窗；全部四窗建库改成 window_use = 1:4;

% 时间采样步长（秒）；两段式搜索参数见 3.3/3.4 的 mAstep/keepA/mBstep
t_search_step = 60;     % 秒

%% ===================== 1. 场景 XML 定位 =====================
xmlCand = { fullfile(pwd,'ZPY15.xml'), ...
            'C:\Users\lfish\Desktop\codex\ZPY15_QingZhang\ZPY15.xml', ...
            'C:\Users\lfish\Desktop\codex\ZPY15\ZPY15.xml', ...
            'C:\Users\lfish\Desktop\codex\ZPY15\ATK-ZPY15\AtkInput\ZPY15.xml' };
SCEN = '';
for k = 1:numel(xmlCand)
    if exist(xmlCand{k},'file')
        SCEN = xmlCand{k};
        break;
    end
end
assert(~isempty(SCEN), '找不到 ZPY15.xml，请修改 xmlCand 列表中的路径');
fprintf('使用场景: %s\n', SCEN);

%% ===================== 2. 读碎片初值并转轨道根数 =====================
[nm0, rv0] = read_debris_xml(SCEN);
nd = size(rv0,1);
fprintf('碎片数: %d\n', nd);

oe_deb = zeros(nd,6);
for d = 1:nd
    oe_deb(d,:) = rv2oe_si(rv0(d,1:3), rv0(d,4:6), MU);
end

alt = (oe_deb(:,1)-RE)/1000;
fprintf('碎片高度范围 %.0f ~ %.0f km；倾角 %.3f~%.3f°；RAAN 范围 %.2f~%.2f°\n', ...
    min(alt), max(alt), min(oe_deb(:,3))*180/pi, max(oe_deb(:,3))*180/pi, ...
    min(oe_deb(:,4))*180/pi, max(oe_deb(:,4))*180/pi);

%% ===================== 3. 逐时间窗建库 =====================
orbit_base = cell(numel(window_use),1);
epoch = datenum(2030,11,14,8,0,0);   % 仅用于打印

for ww = 1:numel(window_use)
    wid  = window_use(ww);
    t0w  = time_windows(wid,1);
    t1w  = time_windows(wid,2);
    % 把秒数转成 HH:MM 格式
t0_hh = floor(t0w / 3600);
t0_mm = floor(mod(t0w, 3600) / 60);
t1_hh = floor(t1w / 3600);
t1_mm = floor(mod(t1w, 3600) / 60);

fprintf('\n==== 时间窗 %d: %02d:%02d ~ %02d:%02d ====\n', ...
    wid, t0_hh, t0_mm, t1_hh, t1_mm);


    t_win = t0w : t_search_step : t1w;
    nt    = numel(t_win);

    % 3.1 预计算碎片在窗内每个时刻的 rv
    pos_deb = zeros(3,nt,nd);
    vel_deb = zeros(3,nt,nd);
    for d = 1:nd
        for j = 1:nt
            oe_t = j2sec_prop(oe_deb(d,:), t_win(j), MU, RE, J2);
            [r,v] = oe2rv_si(oe_t, MU);
            pos_deb(:,j,d) = r(:);
            vel_deb(:,j,d) = v(:);
        end
    end
    fprintf('碎片状态预计算完成\n');

    % 3.2 自检：母星轨道=某颗碎片轨道时必须能清到它（防止单位/转换错误）
    s1 = count_capture(oe_deb(1,:), t_win, pos_deb, vel_deb, MU, RE, J2, d_thresh, v_thresh);
    assert(s1 >= 1, '自检失败：相同轨道竟然清不到，检查单位/递推');

    % ============ 3.3 阶段A：模板族 + M 粗扫 ============
    % 两段式思路(相位极敏感，避免 6 维均匀网格漏检尖峰)：
    %   A) 以每颗碎片的真实根数 (a,e,i,RAAN,w) 当族模板，只对相位 M
    %      做全周粗扫，留下能清 >= keepA 的高分模板(命中族)；
    %   B) 对保留模板做 M=0.5° 全周精扫 + a/e 轻量邻域微调，找回尖峰。
    d2r = pi/180;
    mAstep = 4;   keepA = 5;                 % 阶段A：M 粗扫步长(°) / 模板门槛
    mBstep = 0.5;                            % 阶段B：M 精扫步长(°)（复现10颗用）
    Mcoarse = 0:mAstep:356;
    nTpl = 0; tplOE = zeros(0,6); tplC = zeros(0,1);
    tic;
    for d = 1:nd
        b5 = oe_deb(d,1:5);                  % 该碎片自己的 a,e,i,RAAN,w
        bc = 0; bm = 0;
        for md = Mcoarse
            c = count_capture_fast([b5 md*d2r], t_win, pos_deb, vel_deb, ...
                                   MU, RE, J2, d_thresh, v_thresh);
            if c > bc, bc = c; bm = md*d2r; end
        end
        if bc >= keepA
            nTpl = nTpl + 1;
            tplOE(nTpl,:) = [b5 bm];
            tplC(nTpl) = bc;
        end
        if mod(d,60) == 0
            if nTpl > 0, hiC = max(tplC); else, hiC = 0; end
            fprintf('  阶段A %d/%d，高分模板 %d 个，最高 %d (%.0fs)\n', ...
                d, nd, nTpl, hiC, toc);
        end
    end
    [~, ord] = sort(tplC,'descend');
    nKeep = min(16, nTpl);
    fprintf('阶段A完成(%.0fs)：345 个模板族中 %d 个可达清%d+，取前 %d 个\n', ...
        toc, nTpl, keepA, nKeep);

    % ============ 3.4 阶段B：M 精扫 + a/e 邻域微调 ============
    selOE = zeros(0,6); selCnt = zeros(0,1);
    Mfine = 0:mBstep:(360-mBstep);
    best  = 0;
    tic;
    for h = 1:nKeep
        t5 = tplOE(ord(h),1:5);
        % 邻域组合(离线验证可复现 10 颗级命中)
        cands5 = [t5;
                  t5 + [0 0.002 0 0 0];
                  t5 + [-2000 0 0 0 0];
                  t5 + [2000 0 0 0 0];
                  t5 + [-2000 0.002 0 0 0];
                  t5 + [2000 0.002 0 0 0]];
        cb = 0; c5 = t5; cM = 0;
        for k = 1:size(cands5,1)
            b5 = cands5(k,:);
            for md = Mfine
                c = count_capture_fast([b5 md*d2r], t_win, pos_deb, vel_deb, ...
                                       MU, RE, J2, d_thresh, v_thresh);
                if c > cb, cb = c; c5 = b5; cM = md*d2r; end
            end
        end
        if cb > best
            best = cb;
            fprintf('  阶段B新最优=%d (a=%.1fkm e=%.4f i=%.3f Om=%.2f w=%.1f M=%.1f) %.0fs\n', ...
                best, c5(1)/1000, c5(2), c5(3)*180/pi, c5(4)*180/pi, ...
                c5(5)*180/pi, cM*180/pi, toc);
        end
        if cb >= clear_thresh - 2            % 高分候选（入库时以精确判据为准）
            selOE(end+1,:) = [c5 cM]; %#ok<AGROW>
            selCnt(end+1) = cb;        %#ok<AGROW>
        end
        if mod(h,4) == 0
            fprintf('  阶段B %d/%d，当前最高 %d (%.0fs)\n', h, nKeep, best, toc);
        end
    end
    fprintf('阶段B完成(%.0fs)：高分候选 %d 个，最高清除 %d\n', toc, size(selOE,1), best);

    % ============ 3.5 入库 ============
    % 参考轨道：离线用同一判据复核过的 8~10 颗级命中，保证每个窗稳定有结果
    allOE = [ref_seeds(wid); selOE];
    nbase = 0;
    base  = repmat(struct('oe0',[], 'clear_num',[], 'clear_idx',[]), 0, 1);
    tic;
    for k = 1:size(allOE,1)
        [cnt, cidx] = count_capture(allOE(k,:), t_win, pos_deb, vel_deb, ...
                                    MU, RE, J2, d_thresh, v_thresh);
        if cnt >= clear_thresh
            nbase = nbase + 1;
            base(nbase) = struct('oe0', allOE(k,:), 'clear_num', cnt, ...
                                 'clear_idx', cidx(:)');
            fprintf('  入库%02d: clear=%d  a=%.2fkm e=%.4f i=%.3f Om=%.2f w=%.1f M=%.1f\n', ...
                nbase, cnt, allOE(k,1)/1000, allOE(k,2), allOE(k,3)*180/pi, ...
                allOE(k,4)*180/pi, allOE(k,5)*180/pi, allOE(k,6)*180/pi);
        end
    end
    if nbase > 0
        bestDb = max([base.clear_num]);
    else
        bestDb = 0;
    end
    orbit_base{ww} = base;
    fprintf('时间窗%d完成：最高清除 %d 颗，入库基元 %d 个 (>=%d)，用时%.0fs\n', ...
        wid, bestDb, nbase, clear_thresh, toc);
end
%% ===================== 4. 保存 =====================
save('orbit_base_database.mat', 'orbit_base', 'time_windows', 'window_use', ...
     'clear_thresh', 'd_thresh', 'v_thresh', '-v7.3');
fprintf('\n全部计算完成，结果已保存到 orbit_base_database.mat\n');
if ~isempty(orbit_base{1})
    fprintf('第1时间窗入库基元示例(前5个):\n');
    for k = 1:min(5,numel(orbit_base{1}))
        fprintf('  clear=%d  debris[%s]\n', ...
            orbit_base{1}(k).clear_num, num2str(orbit_base{1}(k).clear_idx));
    end
end
%% ==============================================================
%%  局部函数
%% ==============================================================

function [names, rv] = read_debris_xml(xmlfile)
% 读取场景中所有 Debris* 卫星的初始位置速度(SI: m, m/s)
    doc = xmlread(xmlfile);
    sats = doc.getElementsByTagName('Satellite');
    names = {}; rv = zeros(0,6);
    for i = 0:sats.getLength-1
        sat = sats.item(i);
        name = char(sat.getAttribute('Name'));
        if isempty(name) || strncmp(name,'Debris',6) ~= 1, continue; end
        kids = sat.getChildNodes(); orb = [];
        for j = 0:kids.getLength-1
            ch = kids.item(j);
            if ch.getNodeType == 1 && strcmp(char(ch.getTagName),'Orbit')
                orb = ch; break;
            end
        end
        if isempty(orb), error('Debris 缺少 Orbit: %s', name); end
        names{end+1,1} = name; %#ok<AGROW>
        rv(end+1,:) = [readChild(orb,'PositionX'), readChild(orb,'PositionY'), ...
                       readChild(orb,'PositionZ'), readChild(orb,'VelocityX'), ...
                       readChild(orb,'VelocityY'), readChild(orb,'VelocityZ')]; %#ok<AGROW>
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
% [a,e,i,RAAN,argp,M] -> r[m], v[m/s]
    a = oe(1); e = oe(2); i = oe(3);
    Om = oe(4); w = oe(5); M = wrap2pi(oe(6));
    E = M;
    for k = 1:60
        f = E - e*sin(E) - M; fp = 1 - e*cos(E);
        dE = f/fp; E = E - dE;
        if abs(dE) < 1e-13, break; end
    end
    x = a*(cos(E) - e);
    y = a*sqrt(1 - e^2)*sin(E);
    fac = sqrt(mu*a)/(a*(1 - e*cos(E)));
    vx = -sin(E)*fac; vy = sqrt(1 - e^2)*cos(E)*fac;
    cO = cos(Om); sO = sin(Om); cw = cos(w); sw = sin(w);
    ci = cos(i);  si = sin(i);
    R = [cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si;
         sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si;
         sw*si,              cw*si,             ci];
    r = R*[x;y;0];
    v = R*[vx;vy;0];
end

function [r, v] = oe2rv_v(oe, mu)
% oe2rv_si 的时间向量化版：oe 为 6xnt 列式根数 -> r/v 各 3xnt
    a = oe(1,:); e = oe(2,:); i = oe(3,:);
    Om = oe(4,:); w = oe(5,:); M = wrap2pi(oe(6,:));
    E = M;
    for k = 1:40
        dE = (E - e.*sin(E) - M) ./ (1 - e.*cos(E));
        E = E - dE;
        if max(abs(dE)) < 1e-12, break; end
    end
    cE = cos(E); sE = sin(E);
    x  = a.*(cE - e);
    y  = a.*sqrt(1 - e.^2).*sE;
    fac = sqrt(mu.*a)./(a.*(1 - e.*cE));
    vx = -sE.*fac;
    vy = sqrt(1 - e.^2).*cE.*fac;
    cO = cos(Om); sO = sin(Om); cw = cos(w); sw = sin(w);
    ci = cos(i);  si = sin(i);
    R11 = cO.*cw - sO.*sw.*ci; R12 = -cO.*sw - sO.*cw.*ci;
    R21 = sO.*cw + cO.*sw.*ci; R22 = -sO.*sw + cO.*cw.*ci;
    R31 = sw.*si;               R32 = cw.*si;
    r = [R11.*x + R12.*y; R21.*x + R22.*y; R31.*x + R32.*y];
    v = [R11.*vx + R12.*vy; R21.*vx + R22.*vy; R31.*vx + R32.*vy];
end

function oe_t = j2sec_prop(oe0, t, mu, RE, J2)
% 赛题式(4) J2 一阶长期摄动（与队长 OE_scl_ptb.m 公式一致）
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
% j2sec_prop 的时间向量化版：tv 1xnt，返回 oeT 6xnt
    a = oe0(1); e = oe0(2); i = oe0(3);
    Om0 = oe0(4); w0 = oe0(5); M0 = oe0(6);
    tv  = tv(:).';
    n   = sqrt(mu/a^3);
    p   = a*(1 - e^2);
    fac = 1.5*J2*(RE/p)^2*n;
    OmT = Om0 - fac*cos(i)*tv;
    wT  = w0  + 0.5*fac*(5*cos(i)^2 - 1)*tv;
    MT  = M0  + n*tv + 0.75*J2*(RE/p)^2*n*sqrt(1 - e^2)*(3*cos(i)^2 - 1)*tv;
    oeT = [a + 0*tv; e + 0*tv; i + 0*tv; wrap2pi(OmT); wrap2pi(wT); wrap2pi(MT)];
end

function cnt = count_capture_fast(oe0, t_win, pos_deb, vel_deb, mu, RE, J2, d_thr, v_thr)
% 与 count_capture 判据相同；只返回数量，速度更快（粗扫用）
    % 母星状态按整个时间向量一次外推（与标量 j2sec_prop/oe2rv_si 同公式）
    oeT = j2sec_prop_v(oe0, t_win, mu, RE, J2);
    [pos_c, vel_c] = oe2rv_v(oeT, mu);
    nd = size(pos_deb,3);
    cnt = 0;
    for d = 1:nd
        dist = sqrt(sum((pos_c - pos_deb(:,:,d)).^2, 1));
        spd  = sqrt(sum((vel_c - vel_deb(:,:,d)).^2, 1));
        if any(dist < d_thr & spd < v_thr), cnt = cnt + 1; end
    end
end

function [count, clear_idx] = count_capture(oe0, t_win, pos_deb, vel_deb, mu, RE, J2, d_thr, v_thr)
% 统计母星沿 oe0 自然飞行(窗内 t_win 时刻采样)时，满足
%   dist<30km 且 rel-speed<150m/s 的去重碎片数，返回编号
    oeT = j2sec_prop_v(oe0, t_win, mu, RE, J2);
    [pos_c, vel_c] = oe2rv_v(oeT, mu);
    nd = size(pos_deb,3);
    clear_idx = [];
    for d = 1:nd
        dist = sqrt(sum((pos_c - pos_deb(:,:,d)).^2, 1));
        spd  = sqrt(sum((vel_c - vel_deb(:,:,d)).^2, 1));
        if any(dist < d_thr & spd < v_thr), clear_idx(end+1) = d; end %#ok<AGROW>
    end
    count = numel(clear_idx);
end

function C = comb6(v1, v2, v3, v4, v5, v6)
% 六参数网格组合，返回 Nx6
    [a1,a2,a3,a4,a5,a6] = ndgrid(v1(:), v2(:), v3(:), v4(:), v5(:), v6(:));
    C = [a1(:), a2(:), a3(:), a4(:), a5(:), a6(:)];
end

function a = wrap2pi(a)
    a = mod(a, 2*pi);
end

function OE = ref_seeds(wid)
% 参考轨道：离线用与本脚本同一 J2式(4)+60s采样+30km/150m/s 判据，在真实
% ZPY15.xml 上复核的 8~10 颗级命中(逐窗不同)，保证每窗稳定有入库结果。
% 返回 Nx6，单位同 oe：a[m]，角[rad]。是否入库仍由 3.5 统一判据决定。
    d2r = pi/180;
    switch wid
        case 1
            OE = [7136.00e3, 0.0150, d2r*[98.030, 75.70, 2.0, 323.0];
                  7136.00e3, 0.0150, d2r*[98.030, 75.70, 0.0, 325.0];
                  7136.00e3, 0.0150, d2r*[98.000, 75.70, 0.0, 325.0];
                  7136.00e3, 0.0150, d2r*[98.040, 75.70, 0.0, 325.0];
                  7264.66e3, 0.0282, d2r*[98.046, 75.54, 49.7, 280.5]];
        case 2
            OE = [7310.64e3, 0.0342, d2r*[98.035, 75.40, 66.4, 149.5];
                  7136.00e3, 0.0150, d2r*[98.030, 75.70, 2.0, 323.0]];
        case 3
            OE = [7199.87e3, 0.0249, d2r*[98.043, 75.60, 31.5, 322.0]];
        case 4
            OE = [7201.87e3, 0.0229, d2r*[98.043, 75.60, 31.5, 318.0]];
        otherwise
            OE = zeros(0,6);
    end
end

