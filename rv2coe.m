function [orbital_elements] = rv2coe(r, v, mu)

    if nargin < 3
        mu = 398600.4418;      % 缺省使用地球引力常数
    end

    % 强制转为列向量
    r = r(:); v = v(:);

    r_norm = norm(r);
    v_norm = norm(v);
    v_r = dot(r, v) / r_norm;  % 径向速度（用于真近点角象限判断）

    % 1. 角动量矢量
    h = cross(r, v);
    h_norm = norm(h);

    % 2. 节线矢量 (指向升交点)
    K = [0; 0; 1];
    n_vec = cross(K, h);
    n_norm = norm(n_vec);

    % 3. 偏心率矢量
    e_vec = (1/mu) * ((v_norm^2 - mu/r_norm) * r - dot(r, v) * v);
    e = norm(e_vec);

    % 4. 半长轴
    epsilon = v_norm^2/2 - mu/r_norm;
    a = -mu / (2 * epsilon);

    % 5. 轨道倾角 (用 i 表示)
    i = acos(h(3) / h_norm);

    tol = 1e-12;

    %% 赤道轨道处理
    if abs(h(3)) / h_norm < tol
        Omega = 0;  % 升交点赤经无定义，置 0

        if e < tol  % 圆赤道轨道
            omega = 0;
            M = atan2(r(2), r(1));
            if M < 0, M = M + 2*pi; end
            return;
        else        % 椭圆赤道轨道
            omega = atan2(e_vec(2), e_vec(1));
            % 真近点角 nu
            cos_nu = dot(e_vec, r) / (e * r_norm);
            nu = acos(cos_nu);
            if v_r > 0
                nu = 2*pi - nu;
            end
            % 偏近点角 E -> 平近点角 M
            cos_E = (e + cos_nu) / (1 + e * cos_nu);
            sin_E = (sqrt(1 - e^2) * sin(nu)) / (1 + e * cos_nu);
            E = atan2(sin_E, cos_E);
            M = E - e * sin(E);
            if M < 0, M = M + 2*pi; end
            return;
        end
    end

    %% 一般情况 (非赤道)

    % 6. 升交点赤经 Omega
    Omega = acos(n_vec(1) / n_norm);
    if n_vec(2) < 0
        Omega = 2*pi - Omega;
    end

    %% 圆轨道处理 
    if e < tol
        omega = 0;  % 近地点幅角无定义，置 0
        % 纬度参数 u = omega + nu，因 omega=0，M = u
        u = atan2(r(3) / sin(i), r(1)*cos(Omega) + r(2)*sin(Omega));
        if u < 0, u = u + 2*pi; end
        M = u;
        return;
    end

    % 7. 近地点幅角 omega
    cos_omega = dot(n_vec, e_vec) / (n_norm * e);
    omega = acos(cos_omega);
    if e_vec(3) < 0
        omega = 2*pi - omega;
    end

    % 8. 真近点角 nu
    cos_nu = dot(e_vec, r) / (e * r_norm);
    nu = acos(cos_nu);
    if v_r > 0          % 航天器正在远离近地点
        nu = 2*pi - nu;
    end

    % 9. 偏近点角 E（无奇异公式）
    cos_E = (e + cos_nu) / (1 + e * cos_nu);
    sin_E = (sqrt(1 - e^2) * sin(nu)) / (1 + e * cos_nu);
    E = atan2(sin_E, cos_E);

    % 10. 平近点角 M (开普勒方程)
    M = E - e * sin(E);
    if M < 0
        M = M + 2*pi;
    end

    % 归一化角度到 [0, 2*pi)
    if omega < 0, omega = omega + 2*pi; end
    if Omega < 0, Omega = Omega + 2*pi; end

    orbital_elements = [a, e, i, Omega, omega, M];
end