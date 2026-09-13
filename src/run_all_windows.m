%% run_all_windows.m
clc; clear;

global J_2 R_E_m mu
J_2 = 1.082627e-3;  R_E_m = 6378137;  mu = 3.986004418e14;

S = load('result_24h.mat');
Dstate_t1_ra = S.Dstate_t1_ra;
Dstate_t1_va = S.Dstate_t1_va;

% 与冠军方案一致的 4 个时段
windows = [
    0,  8;    % 08:00-16:00
    6, 14;    % 14:00-22:00
   12, 20;    % 20:00-04:00(+1)
   16, 24;    % 00:00-08:00(+1)
];

all_segments = cell(size(windows, 1), 1);
for w = 1:size(windows, 1)
    tag = sprintf('W%d', w);
    all_segments{w} = search_window(windows(w,1), windows(w,2), ...
        Dstate_t1_ra, Dstate_t1_va, mu, J_2, R_E_m, tag);
end

save('all_segments.mat', 'all_segments', 'windows', '-v7.3');
fprintf('4 个时段数据库已生成\n');
for w = 1:4
    ncap = [all_segments{w}.n_captured];
    fprintf('  W%d (%.0fh-%.0fh)：%d 条，≥8 颗 %d 条\n', ...
        w, windows(w,1), windows(w,2), numel(ncap), sum(ncap >= 8));
end