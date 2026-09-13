%% ===================== 组合优化主程序 =====================
% 数据库字段：oe0（轨道六根数）, clear_num（该轨道清除碎片总数）, clear_idx（清除碎片编号列表）
% 算法结构：局部优化（等优接受）和全局优化（窗口级交叉）交替执行
% 增强：多次重启 + 精英池 + 自适应 p_local
% 依赖：先跑 deb_clear.m 生成 orbit_base_database.mat
% ================================================================================
clear; clc; close all;

% -------- 参数 --------
MU = 3.986004418e14; RE = 6378137; J2 = 1.08262668e-3;
d_thresh = 30000; v_thresh = 150; t_search_step = 60;
transition_time = 45*60;
n_sc = 3;

% -------- 加载数据库 --------
assert(exist('orbit_base_database.mat','file')==2, '请先运行 deb_clear.m');
S = load('orbit_base_database.mat');
orbit_base   = S.orbit_base;
time_windows = S.time_windows;
window_use   = S.window_use;
n_win        = numel(window_use);
for ww = 1:n_win
    fprintf('窗口%d: %d 基元\n', ww, numel(orbit_base{ww}));
    % 校验字段
    assert(isfield(orbit_base{ww}, 'oe0'),        '字段名应为 oe0');
    assert(isfield(orbit_base{ww}, 'clear_num'),  '字段名应为 clear_num');
    assert(isfield(orbit_base{ww}, 'clear_idx'),  '字段名应为 clear_idx');
end

% -------- 读取碎片 --------
assert(exist('ZPY15.xml','file')==2, '未找到 ZPY15.xml');
[~, rv0] = read_debris_xml('ZPY15.xml');
nd = size(rv0,1);
oe_deb = zeros(nd,6);
for d = 1:nd
    oe_deb(d,:) = rv2oe_si(rv0(d,1:3), rv0(d,4:6), MU);
end
fprintf('碎片总数: %d\n', nd);

% -------- 预计算每个基元的清除时刻 --------
precalc = cell(n_win,1);
window_end_t = zeros(n_win,1);
for ww = 1:n_win
    wid = window_use(ww);
    t0w = time_windows(wid,1); t1w = time_windows(wid,2);
    window_end_t(ww) = t1w;
    t_win = t0w:t_search_step:t1w;
    nt = numel(t_win);
    pos_d = zeros(nt,3,nd); vel_d = zeros(nt,3,nd);
    for d = 1:nd
        oeT = j2sec_prop_v(oe_deb(d,:), t_win, MU, RE, J2);
        [r,v] = oe2rv_v(oeT, MU);
        pos_d(:,:,d) = r'; vel_d(:,:,d) = v';
    end
    nprim = numel(orbit_base{ww});
    pc = cell(nprim,1);
    for i = 1:nprim
       [ids, times] = find_capture_times(orbit_base{ww}(i).oe0, ...
            t_win, pos_d, vel_d, MU, RE, J2, d_thresh, v_thresh);
        pc{i}.ids = ids; pc{i}.times = times;
    end
    precalc{ww} = pc;
end

obj_fun = @(sol) eval_solution(sol, precalc, n_sc, n_win, ...
                                window_end_t, transition_time, nd);

%% ======================= 多次重启 + 精英池 =======================
n_restarts     = 8;
pop_size       = 30;
max_iter       = 8000;
max_no_improve = 1500;
elite_size     = 5;

best_global_obj = -inf;
best_global_sol = [];
elite_pool      = {};
history_all     = [];

for restart = 1:n_restarts
    rng(1000 + restart);
    fprintf('\n######## 重启 %d/%d ########\n', restart, n_restarts);

    % ---------- 初始化种群（精英池注入） ----------
    population = cell(pop_size,1);
    pop_obj    = zeros(pop_size,1);
    n_elite_use = min(elite_size, numel(elite_pool));
    for i = 1:n_elite_use
        population{i} = elite_pool{i};
        pop_obj(i)    = obj_fun(population{i});
    end
    for i = n_elite_use+1:pop_size
        population{i} = random_solution(orbit_base, n_sc, n_win);
        pop_obj(i)    = obj_fun(population{i});
    end
    pop_obj = pop_obj(:);

    best_obj   = max(pop_obj);
    no_improve = 0;
    history    = zeros(max_iter, 1);
    fprintf('初始: best=%d\n', best_obj);

    % ---------- 主循环：交替执行 ----------
    for iter = 1:max_iter
        p_local = 0.4 + 0.3 * (iter / max_iter);   % 自适应

    if rand() < p_local
            % ============ 局部优化：随机换一个段，等优即接受 ============
            i = randi(pop_size);            % 随机选一个个体
            sol = population{i};            % 把这个解放进 sol
            obj = pop_obj(i);               % 记下它当前的分数

            k = randi(n_sc);                % 随机选一个航天器 (1/2/3)
            w = randi(n_win);               % 随机选一个窗口   (1/2/3/4)

            nprim = numel(orbit_base{w});   % 看这个窗口总共有几个基元
            cands = setdiff(1:nprim, sol(:,w));   %  找出这个窗口"还没被用"的基元

        if ~isempty(cands)              % 如果还有没用过的
        sol(k,w) = cands(randi(numel(cands)));   % 随机换一个
        new_obj  = obj_fun(sol);                 % 算新分数
        if new_obj >= obj                        % 新分 >= 旧分？
        population{i} = sol;                 % 是 → 接受
        pop_obj(i)    = new_obj;
        end                                      % 否 → 什么都不做（自动丢弃）
    end
        else
            % ============ 全局优化：窗口级交叉 ============
            i1 = randi(pop_size);                  %  选第一个父代
            i2 = randi(pop_size);                  %  选第二个父代
        while i2 == i1                         % 确保它们不同
            i2 = randi(pop_size);
        end

    [c1, c2] = crossover_window(...);      % 交叉，生成两个子代
    o1 = obj_fun(c1);                      % 算子代1的分数
    o2 = obj_fun(c2);                      % 算子代2的分数

    [~, wi] = min(pop_obj);                % 找当前种群最差的个体
    if o1 > pop_obj(wi)                    % 子代1比最差的强？
    population{wi}=c1; pop_obj(wi)=o1; % 替换掉最差的
    end

    [~, wi] = min(pop_obj);                % 再找当前最差的
    if o2 > pop_obj(wi)                    % 子代2比最差的强？
    population{wi}=c2; pop_obj(wi)=o2; % ⑫ 替换掉
    end

        % ---------- 收敛判据 ----------
        cur_best = max(pop_obj);
        if cur_best > best_obj
            best_obj   = cur_best;
            no_improve = 0;
        else
            no_improve = no_improve + 1;
        end
        history(iter) = best_obj;

        if no_improve >= max_no_improve
            fprintf('  迭代 %d 早停\n', iter);
            history = history(1:iter);
            break;
        end
    end
    history_all = [history_all; history]; 

    % ---------- 更新精英池 ----------
    [~, ord] = sort(pop_obj, 'descend');
    for ii = 1:min(elite_size, numel(ord))
        sol_i = population{ord(ii)};
        is_dup = false;
        for e = 1:numel(elite_pool)
            if isequal(elite_pool{e}, sol_i), is_dup=true; break; end
        end
        if ~is_dup
            elite_pool{end+1} = sol_i; 
        end
    end
    if numel(elite_pool) > elite_size
        ep_obj = zeros(numel(elite_pool),1);
        for e = 1:numel(elite_pool)
            ep_obj(e) = obj_fun(elite_pool{e});
        end
        [~, eo] = sort(ep_obj, 'descend');
        elite_pool = elite_pool(eo(1:elite_size));
    end

    fprintf('  重启%d结束: best=%d (iter=%d)\n', restart, best_obj, iter);

    if best_obj > best_global_obj
        best_global_obj = best_obj;
        [~, bi] = max(pop_obj);
        best_global_sol = population{bi};
        fprintf('  *** 全局最优更新: %d ***\n', best_global_obj);
    end
end

%% ======================= 输出结果 =======================
best_sol   = best_global_sol;
final_best = best_global_obj;

fprintf('\n================= 全局最优方案 =================\n');
fprintf('清除碎片总数: %d\n\n', final_best);

all_deb_global = false(1, nd);
for k = 1:n_sc
    sc_deb = false(1, nd);
    fprintf('航天器 %d:\n', k);
    for ww = 1:n_win
        idx   = best_sol(k, ww);
        p     = orbit_base{ww}(idx);
        ids   = precalc{ww}{idx}.ids;
        times = precalc{ww}{idx}.times;
        if ww < n_win
            tmax = window_end_t(ww) - transition_time;
            ids  = ids(times <= tmax);
        end
        sc_deb(ids) = true;
       
        fprintf('  窗%d: 基元#%-4d N_p=%2d  a=%.1f km  e=%.4f  i=%.3f°  M=%.1f°\n', ...
            ww, idx, p.clear_num, p.oe0(1)/1000, p.oe0(2), ...
            p.oe0(3)*180/pi, p.oe0(6)*180/pi);
    end
    fprintf('  → 有效清除: %d\n\n', sum(sc_deb));
    all_deb_global = all_deb_global | sc_deb;
end
fprintf('全局去重后清除: %d 个碎片\n', sum(all_deb_global));

save('combine_result_v5.mat', 'best_sol', 'final_best', ...
     'elite_pool', 'history_all', '-v7.3');
fprintf('已保存 combine_result_v5.mat\n');

% ---------- 收敛曲线 ----------
figure('Color','w');
plot(history_all, 'b-', 'LineWidth', 1.2);
xlabel('迭代次数（跨重启累计）'); ylabel('清除碎片数');
title('组合优化收敛曲线'); grid on;

%% ===================== 辅助函数 =====================

% ---------- 随机解 ----------
function sol = random_solution(orbit_base, n_sc, n_win)
    sol = zeros(n_sc, n_win);
    for ww = 1:n_win
        nprim = numel(orbit_base{ww});
        sol(:,ww) = randperm(nprim, n_sc)';
    end
end

% ---------- 目标函数 ----------
function score = eval_solution(sol, precalc, n_sc, n_win, ...
                                window_end_t, transition_time, nd)
    all_deb = false(1, nd);
    for k = 1:n_sc
        sc_deb = false(1, nd);
        for ww = 1:n_win
            idx   = sol(k, ww);
            ids   = precalc{ww}{idx}.ids;
            times = precalc{ww}{idx}.times;
            if ww < n_win
                tmax = window_end_t(ww) - transition_time;
                ids  = ids(times <= tmax);
            end
            if ~isempty(ids), sc_deb(ids) = true; end
        end
        all_deb = all_deb | sc_deb;
    end
    score = sum(all_deb);
end

% ---------- 窗口级交叉：整列交换 ----------
function [c1, c2] = crossover_window(p1, p2, n_win)
    c1 = p1;                           % 子代1先复制父代1
    c2 = p2;                           % 子代2先复制父代2
    if n_win > 1
        cut = randi(n_win - 1);        % 随机选切点（1、2 或 3）
        c1(:, cut+1:end) = p2(:, cut+1:end);  % 子代1的后半段来自父代2
        c2(:, cut+1:end) = p1(:, cut+1:end);  % 子代2的后半段来自父代1
    end
end

% ---------- 求某基元在窗口内首次清除每个碎片的时刻 ----------
function [ids, times] = find_capture_times(oe0, t_win, pos_d, vel_d, ...
                                            mu, RE, J2, d_thr, v_thr)
    nt = numel(t_win); nd = size(pos_d,3);
    oeT = j2sec_prop_v(oe0, t_win, mu, RE, J2);
    [pc, vc] = oe2rv_v(oeT, mu);
    pc = pc'; vc = vc';
    dsq = sum((pc - pos_d).^2, 2);
    vsq = sum((vc - vel_d).^2, 2);
    hit = reshape(dsq < d_thr^2 & vsq < v_thr^2, nt, nd);
    ids = zeros(0,1); times = zeros(0,1);
    for d = 1:nd
        ti = find(hit(:,d), 1, 'first');
        if ~isempty(ti)
            ids(end+1,1)   = d;          %#ok<AGROW>
            times(end+1,1) = t_win(ti);  %#ok<AGROW>
        end
    end
end

% ---------- 以下函数与 deb_clear.m 一致 ----------
function [names, rv] = read_debris_xml(xmlfile)
    doc  = xmlread(xmlfile);
    sats = doc.getElementsByTagName('Satellite');
    names = {}; rv = zeros(0,6);
    for i = 0:sats.getLength-1
        sat  = sats.item(i);
        name = char(sat.getAttribute('Name'));
        if isempty(name) || ~startsWith(name,'Debris'), continue; end
        kids = sat.getChildNodes(); orb = [];
        for j = 0:kids.getLength-1
            ch = kids.item(j);
            if ch.getNodeType==1 && strcmp(char(ch.getTagName),'Orbit')
                orb = ch; break;
            end
        end
        if isempty(orb), error('Debris 缺少 Orbit: %s', name); end
        names{end+1,1} = name;  %#ok<AGROW>
        rv(end+1,:) = [readChild(orb,'PositionX'), readChild(orb,'PositionY'), ...
                       readChild(orb,'PositionZ'), readChild(orb,'VelocityX'), ...
                       readChild(orb,'VelocityY'), readChild(orb,'VelocityZ')];  %#ok<AGROW>
    end
end
function val = readChild(node, tag)
    kids = node.getChildNodes(); val = NaN;
    for j = 0:kids.getLength-1
        ch = kids.item(j);
        if ch.getNodeType==1 && strcmp(char(ch.getTagName),tag)
            val = str2double(char(ch.getTextContent)); return;
        end
    end
end
function oe = rv2oe_si(r, v, mu)
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
        u = atan2(r(3)/sin(inc), r(1)*cos(RAAN)+r(2)*sin(RAAN));
        M = wrap2pi(u);
    else
        argp = acos(max(-1,min(1,dot(nv,evec)/(nn*e))));
        if evec(3) < 0, argp = 2*pi - argp; end
        cnu = max(-1,min(1,dot(evec,r)/(e*rn)));
        nu  = acos(cnu);
        if dot(r,v) < 0, nu = 2*pi - nu; end
        E = 2*atan2(sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2));
        M = wrap2pi(E - e*sin(E));
    end
    oe = [a, e, inc, wrap2pi(RAAN), wrap2pi(argp), M];
end
function [r, v] = oe2rv_v(oe, mu)
    a = oe(1,:); e = oe(2,:); i = oe(3,:);
    Om = oe(4,:); w = oe(5,:); M = wrap2pi(oe(6,:));
    E = M;
    for k = 1:40
        dE = (E - e.*sin(E) - M) ./ (1 - e.*cos(E));
        E = E - dE;
        if max(abs(dE)) < 1e-12, break; end
    end
    cE = cos(E); sE = sin(E);
    x = a.*(cE - e); y = a.*sqrt(1 - e.^2).*sE;
    fac = sqrt(mu.*a)./(a.*(1 - e.*cE));
    vx = -sE.*fac; vy = sqrt(1 - e.^2).*cE.*fac;
    cO = cos(Om); sO = sin(Om); cw = cos(w); sw = sin(w);
    ci = cos(i); si = sin(i);
    R11 = cO.*cw - sO.*sw.*ci; R12 = -cO.*sw - sO.*cw.*ci;
    R21 = sO.*cw + cO.*sw.*ci; R22 = -sO.*sw + cO.*cw.*ci;
    R31 = sw.*si; R32 = cw.*si;
    r = [R11.*x + R12.*y; R21.*x + R22.*y; R31.*x + R32.*y];
    v = [R11.*vx + R12.*vy; R21.*vx + R22.*vy; R31.*vx + R32.*vy];
end
function oeT = j2sec_prop_v(oe0, tv, mu, RE, J2)
    a = oe0(1); e = oe0(2); i = oe0(3);
    Om0 = oe0(4); w0 = oe0(5); M0 = oe0(6);
    tv = tv(:).';
    n = sqrt(mu/a^3);
    p = a*(1 - e^2);
    fac = 1.5*J2*(RE/p)^2*n;
    OmT = Om0 - fac*cos(i)*tv;
    wT  = w0 + 0.5*fac*(5*cos(i)^2 - 1)*tv;
    MT  = M0 + n*tv + 0.75*J2*(RE/p)^2*n*sqrt(1-e^2)*(3*cos(i)^2-1)*tv;
    oeT = [a + 0*tv; e + 0*tv; i + 0*tv; wrap2pi(OmT); wrap2pi(wT); wrap2pi(MT)];
end
function a = wrap2pi(a)
    a = mod(a, 2*pi);
end
