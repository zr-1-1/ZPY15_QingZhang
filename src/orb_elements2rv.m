function state = orb_elements2rv(elements, mu)
% ORB_ELEMENTS2RV  由轨道六根数计算位置和速度矢量（ECI坐标系）
%   输入：
%       elements - Nx6 矩阵，每行为 [a, e, i, Omega, omega, M]
%                  a   - 半长轴 [km]
%                  e   - 离心率 [无量纲]
%                  i   - 轨道倾角 [rad]
%                  Omega - 升交点赤经 [rad]
%                  omega - 近地点幅角 [rad]
%                  M   - 平近点角 [rad]
%       mu       - 引力常数 [km^3/s^2]（可选，默认地球常数 398600.4418 km^3/s^2）
%   输出：
%       state    - Nx6 矩阵，每行为 [x, y, z, vx, vy, vz]
%                  位置单位为 km，速度单位为 km/s（ECI惯性系）
%
%   示例：
%       elem = [42164, 0.001, 0.1, 0.5, 0.3, 1.2];   % 注意 a 已是 km
%       state = orb_elements2rv(elem);
%       fprintf('位置 (km): [%.3f, %.3f, %.3f]\n', state(1:3));
%       fprintf('速度 (km/s): [%.6f, %.6f, %.6f]\n', state(4:6));

    % 导入全局变量，保持参数一致性
    global J_2 R_E R_E_m w_E mu

    % 这里没有nargin这个变量吧？为什么要写这个？
    % if nargin < 2
    %     mu = 398600.4418;   % 地球引力常数 (km^3/s^2)
    % end

    % 检查输入是否为 6 列
    if size(elements, 2) ~= 6
        error('elements 必须为 6 列矩阵');
    end

    N = size(elements, 1);   % 获取轨道组数
    a = elements(:,1);
    e = elements(:,2);
    i = elements(:,3);
    Omega = elements(:,4);
    omega = elements(:,5);
    M = elements(:,6);

    state = zeros(N, 6);

    for k = 1:N
        % 1. 解开普勒方程 M = E - e*sin(E)
        if e(k) < 1e-12
            E = M(k);
        else
            E0 = M(k);
            tol = 1e-12;
            maxIter = 100;
            for iter = 1:maxIter
                Ek = E0 - (E0 - e(k)*sin(E0) - M(k)) / (1 - e(k)*cos(E0));
                if abs(Ek - E0) < tol
                    break;
                end
                E0 = Ek;
            end
            E = Ek;
        end

        % 2. 计算真近点角 f 和径向距离 r
        cosE = cos(E);
        sinE = sin(E);
        r = a(k) * (1 - e(k)*cosE);          % 单位 km
        cosf = (cosE - e(k)) / (1 - e(k)*cosE);
        sinf = sqrt(1 - e(k)^2) * sinE / (1 - e(k)*cosE);

        % 3. 轨道平面 (PQW) 中的位置和速度
        p = a(k) * (1 - e(k)^2);              % 单位 km
        r_PQW = [r*cosf; r*sinf; 0];
        v_PQW = sqrt(mu/p) * [-sinf; e(k)+cosf; 0];   % 单位 km/s

        % 4. 旋转至 ECI 坐标系
        cosO = cos(Omega(k)); sinO = sin(Omega(k));
        cosi = cos(i(k));    sini = sin(i(k));
        cosw = cos(omega(k)); sinw = sin(omega(k));

        R = [cosO*cosw - sinO*sinw*cosi, -cosO*sinw - sinO*cosw*cosi, sinO*sini;
             sinO*cosw + cosO*sinw*cosi, -sinO*sinw + cosO*cosw*cosi, -cosO*sini;
             sinw*sini, cosw*sini, cosi];

        r_vec = R * r_PQW;
        v_vec = R * v_PQW;

        state(k,:) = [r_vec', v_vec'];
    end
end