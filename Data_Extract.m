clc;
clear;
J_2 = 1.0827e-3;
R_E = 6378.14;
%% 
conID = atkOpen();
% 获取所有碎片初始轨道根数
Debris_oe = zeros(345,6);
for i = 1:345
    paramstr = sprintf('*/Satellite/Debris%d "14 Nov 2030 08:00:00.000"', i);
    Debris_d = atkConnect(conID, 'Position', paramstr);
    Debris_rv = str2double(split(Debris_d))/1000;
    Debris_oe(i,:) = rv2coe(Debris_rv(1:3),Debris_rv(4:6));
end
atkclose(conID);

%% 按照一定步长生成24小时内碎片每个时刻的轨道根数
t_step = 1;
outercell = cell(1,345);
for i = 1:345
    innercell = cell(1,20000);
    E_0 = Debris_oe(i, :);
    idx = 0;

    for t = 0: t_step: 60*60*24
        E_t = OE_scl_ptb(E_0, t);
        idx = idx+1;
        innercell{idx} = E_t;
    end
    D_E_i = cat(2,innercell{:});
    outercell{i} = D_E_i;
end
D_Eall = cat(3,outercell{:});

