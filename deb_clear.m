%% ===================== 主程序：轨道基元数据库生成（论文 2.1 节） =====================
clear; clc; close all;

%% -------- 参数 --------
MU = 3.986004418e14; RE = 6378137; J2 = 1.08262668e-3;
d_thresh = 30000; v_thresh = 150;
transition_time = 45*60;
n_sc = 3;

time_windows = [0,       8*3600;
                6*3600,  14*3600;
                12*3600, 20*3600;
                16*3600, 24*3600];

%% -------- 1. 加载数据库 --------
assert(exist('high_capture_segments.mat','file')==2, '未找到数据库');
S = load('high_capture_segments.mat');
high_all = S.high_all;
n_win = numel(high_all);
for ww = 1:n_win
    fprintf('窗口%d: %d 基元\n', ww, numel(high_all{ww}));
end
nd = 345;
window_start_t = time_windows(1:n_win, 1);
window_end_t   = time_windows(1:n_win, 2);

%% -------- 2. 预计算 --------
precalc = cell(n_win, 1);
for ww = 1:n_win
    nprim = numel(high_all{ww});
    pc = cell(nprim, 1);
    for i = 1:nprim
        p = high_all{ww}(i);
        ids   = p.captured_ids(:);
        times = p.capture_times(:);
        valid = times > 0;
        pc{i}.ids   = ids(valid);
        pc{i}.times = times(valid);
    end
    precalc{ww} = pc;
end

%% -------- 3. 目标函数 --------
obj_fun = @(sol) eval_chain(sol, precalc, n_sc, n_win, ...
                             window_start_t, window_end_t, ...
                             transition_time, nd);

%% -------- 4. 遗传算法 --------
n_restarts     = 8;
pop_size       = 30;
max_iter       = 8000;
max_no_improve = 1500;
elite_size     = 5;

best_global_obj = -inf;
best_global_sol = [];
elite_pool      = {};

for restart = 1:n_restarts
    rng(1000 + restart);
    fprintf('\n######## 重启 %d/%d ########\n', restart, n_restarts);

    population = cell(pop_size, 1);
    pop_obj    = zeros(pop_size, 1);
    n_elite_use = min(elite_size, numel(elite_pool));
    for i = 1:n_elite_use
        population{i} = elite_pool{i};
        pop_obj(i)    = obj_fun(population{i});
    end
    for i = n_elite_use+1:pop_size
        population{i} = random_solution(high_all, n_sc, n_win);
        pop_obj(i)    = obj_fun(population{i});
    end
    pop_obj = pop_obj(:);

    best_obj   = max(pop_obj);
    no_improve = 0;
    fprintf('初始: best=%d\n', best_obj);

    for iter = 1:max_iter
        p_local = 0.4 + 0.3 * (iter / max_iter);

        if rand() < p_local
            i = randi(pop_size);
            sol = population{i};
            obj = pop_obj(i);
            k = randi(n_sc);
            w = randi(n_win);
            nprim = numel(high_all{w});
            cands = setdiff(1:nprim, sol(:,w));
            if ~isempty(cands)
                sol(k,w) = cands(randi(numel(cands)));
                new_obj  = obj_fun(sol);
                if new_obj >= obj
                    population{i} = sol;
                    pop_obj(i)    = new_obj;
                end
            end
        else
            i1 = randi(pop_size);
            i2 = randi(pop_size);
            while i2 == i1, i2 = randi(pop_size); end
            [c1, c2] = crossover_window(population{i1}, population{i2}, n_win);
            o1 = obj_fun(c1); o2 = obj_fun(c2);
            [~, wi] = min(pop_obj);
            if o1 > pop_obj(wi), population{wi}=c1; pop_obj(wi)=o1; end
            [~, wi] = min(pop_obj);
            if o2 > pop_obj(wi), population{wi}=c2; pop_obj(wi)=o2; end
        end

        cur_best = max(pop_obj);
        if cur_best > best_obj
            best_obj   = cur_best;
            no_improve = 0;
        else
            no_improve = no_improve + 1;
        end
        if no_improve >= max_no_improve
            break;
        end
    end

    [~, ord] = sort(pop_obj, 'descend');
    for ii = 1:min(elite_size, numel(ord))
        sol_i = population{ord(ii)};
        is_dup = false;
        for e = 1:numel(elite_pool)
            if isequal(elite_pool{e}, sol_i), is_dup=true; break; end
        end
        if ~is_dup, elite_pool{end+1} = sol_i; end %#ok<AGROW>
    end
    if numel(elite_pool) > elite_size
        ep_obj = zeros(numel(elite_pool),1);
        for e = 1:numel(elite_pool)
            ep_obj(e) = obj_fun(elite_pool{e});
        end
        [~, eo] = sort(ep_obj, 'descend');
        elite_pool = elite_pool(eo(1:elite_size));
    end

    fprintf('  重启%d: best=%d\n', restart, best_obj);

    if best_obj > best_global_obj
        best_global_obj = best_obj;
        [~, bi] = max(pop_obj);
        best_global_sol = population{bi};
        fprintf('  *** 全局最优: %d ***\n', best_global_obj);
    end
end

%% -------- 5. 输出 --------
best_sol = best_global_sol;

fprintf('\n================= 最优方案 =================\n');

all_deb_global = false(1, nd);
for k = 1:n_sc
    sc_deb = false(1, nd);
    fprintf('\n航天器 %d:\n', k);
    t_prev_end = -inf;

    for ww = 1:n_win
        idx = best_sol(k, ww);
        p   = high_all{ww}(idx);
        ids   = precalc{ww}{idx}.ids;
        times = precalc{ww}{idx}.times;

        times_abs = times;   % 绝对时间

        if ww == 1
            t_start_abs = window_start_t(ww);
        else
            t_start_abs = max(t_prev_end + transition_time, window_start_t(ww));
        end

        keep = (times_abs >= t_start_abs) & (times_abs <= window_end_t(ww));
        ids_v   = ids(keep);
        times_v = times_abs(keep);

        n_orig  = length(ids);
        n_valid = length(ids_v);
        n_lost  = n_orig - n_valid;

        if ~isempty(ids_v)
            sc_deb(ids_v) = true;
            t_prev_end = max(times_v);
        else
            t_prev_end = t_start_abs;
        end

        fprintf('  段%d: 基元#%-4d  原始%2d  → 有效%2d (扣%d)\n', ...
            ww, idx, n_orig, n_valid, n_lost);
        fprintf('        段起点=%.2fh, 段终点=%.2fh\n', ...
            t_start_abs/3600, window_end_t(ww)/3600);
        if ~isempty(times_v)
            fprintf('        捕获区间=%.2fh~%.2fh\n', ...
                min(times_v)/3600, max(times_v)/3600);
        end
    end
    fprintf('  → 有效清除: %d\n', sum(sc_deb));
    all_deb_global = all_deb_global | sc_deb;
end

fprintf('\n==========================================\n');
fprintf('清除碎片总数: %d\n', sum(all_deb_global));
fprintf('==========================================\n');

save('combine_result_chain.mat', 'best_sol', 'all_deb_global', ...
     'best_global_obj', '-v7.3');
fprintf('已保存 combine_result_chain.mat\n');

%% ===================== 辅助函数 =====================

function score = eval_chain(sol, precalc, n_sc, n_win, ...
                             window_start_t, window_end_t, ...
                             transition_time, nd)
    all_deb = false(1, nd);

    for k = 1:n_sc
        sc_deb = false(1, nd);
        t_prev_end = -inf;

        for ww = 1:n_win
            idx   = sol(k, ww);
            ids   = precalc{ww}{idx}.ids;
            times = precalc{ww}{idx}.times;

            times_abs = times;   

            if ww == 1
                t_start_abs = window_start_t(ww);
            else
                t_start_abs = max(t_prev_end + transition_time, ...
                                  window_start_t(ww));
            end

            keep = (times_abs >= t_start_abs) & (times_abs <= window_end_t(ww));
            ids_v   = ids(keep);
            times_v = times_abs(keep);

            if ~isempty(ids_v)
                sc_deb(ids_v) = true;
                t_prev_end = max(times_v);
            else
                t_prev_end = t_start_abs;
            end
        end

        all_deb = all_deb | sc_deb;
    end
    score = sum(all_deb);
end

function sol = random_solution(high_all, n_sc, n_win)
    sol = zeros(n_sc, n_win);
    for ww = 1:n_win
        nprim = numel(high_all{ww});
        sol(:,ww) = randperm(nprim, n_sc)';
    end
end

function [c1, c2] = crossover_window(p1, p2, n_win)
    c1 = p1; c2 = p2;
    if n_win > 1
        cut = randi(n_win - 1);
        c1(:, cut+1:end) = p2(:, cut+1:end);
        c2(:, cut+1:end) = p1(:, cut+1:end);
    end
end
