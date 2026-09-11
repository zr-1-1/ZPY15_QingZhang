%%
clc; clear;
mu = 3.986004418e14; % 地球引力常数，单位：(m^3)/(s^2)
t0 = 0; t1 = 28800; t2 = 57600; t3 = 86400;
t_step = 3; % 时间步长 单位：秒
%% 按照一定步长生成24小时内碎片每个时刻的轨道根数

outercell_r = cell(1,345);
outercell_v = cell(1,345);
for i = 1:345
    innercell_r = cell(1,10000);
    innercell_v = cell(1,10000);
    E_0 = Debris_oe(i, :);
    idx = 0;

    for t = t0: t_step: t1
        E_t = OE_scl_ptb(E_0, t);
        [Dbstate_t1_r, Dbstate_t1_v] = orb_elements2rv(E_t, mu);
        idx = idx+1;
        innercell_r{idx} = Dbstate_t1_r;
        innercell_v{idx} = Dbstate_t1_v;
    end
    
    Dstate_t1_r = horzcat(innercell_r{1:idx});
    Dstate_t1_v = horzcat(innercell_v{1:idx});
    outercell_r{i} = Dstate_t1_r;
    outercell_v{i} = Dstate_t1_v;
end
Dstate_t1_ra = cat(3,outercell_r{:});
Dstate_t1_va = cat(3,outercell_v{:});

clear innercell_v innercell_r outercell_v outercell_r Dstate_t1_r Dstate_t1_v
save("result.mat","Dstate_t1_va","Dstate_t1_ra");

%% 遍历搜索符合条件的轨道段

% 固定参数网格 (不用循环的维度)
a_list = [7000000, 7250000, 7500000, 7750000, 8000000];
i_list = [0, 0.05, 0.10, 0.16];
Omega_list = [1.29, 1.33, 1.37];
e_list = [0, 0.03, 0.06, 0.10, 0.15];
% omega 特殊处理，分段取值
omega_list = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0, 5.5, 6.0, 6.2832];
M_list = 0:0.3:(2*pi);  % 21个
 
% 预分配内存
orb_cell = cell(10000,1);
total_idx = 0;

for idx_Omega = 1:length(Omega_list)
    Omega = Omega_list(idx_Omega);
     for idx_a = 1:length(a_list)
        a = a_list(idx_a);
        for idx_i = 1:length(i_list)
            i = i_list(idx_i);
            for idx_e = 1:length(e_list)
                e = e_list(idx_e);
                for idx_om = 1:length(omega_list)
                    omega = omega_list(idx_om);
                    for idx_M = 1:length(M_list)
                        M = M_list(idx_M);
                        
                        oe = [a, e, i, Omega, omega, M];
                        [Sat_r,Sat_v] = orb_elements2rv(oe,mu);
                        state = [Sat_r; Sat_v];
                        
                        %  调用 J2 积分器 (8小时)
                        duation = t0:t_step:t1; % 和碎片的时刻保持一致
                        opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
                        [~, Y] = ode78(@Sat_J2ODE, duation, state, opts);
                        sc_pos = Y(:, 1:3);
                        sc_vol = Y(:, 4:6);
                
                        pos_re = squeeze(sqrt(sum((sc_pos'-Dstate_t1_ra).^2,1))); % 相对距离
                        vol_re = squeeze(sqrt(sum((sc_vol'-Dstate_t1_va).^2,1))); % 相对速度
                        count = 0;
                        final_t =0;
                        tem = zeros(1,16);
                        for i = 1:1:length(pos_re)
                            for j = 1:1:345
                                if pos_re(i,j)<30000 && vol_re(i,j)<150
                                    count = count+1; % 捕获数
                                    tem(8+count) = j;   % 记录哪个碎片被捕获
                                    final_t = i;     % 最后一个碎片的捕获时刻
                                end
                            end
                        end
             
                            tem(1:6) = oe;
                            tem(7) = count; tem(8) = final_t; % 数组第7个数是捕获数，第8个数是最后一个碎片捕获时间
                            total_idx = total_idx+1;
                            orb_cell{total_idx} = tem;
                        
                    end
                end
            end
        end
    end
end