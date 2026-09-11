function [r, v] = orb_elements2rv(oe, mu)
% orbital_elements_to_state  根据轨道六根数计算位置和速度矢量

a = oe(1); e = oe(2); i = oe(3); RAAN = oe(4); omega = oe(5); M = oe(6);

% 1. 解开普勒方程: M = E - e*sin(E)  (牛顿法)
E = M;  % 初始猜测
tol = 1e-12;
maxIter = 100;

for iter = 1:maxIter
    f = E - e*sin(E) - M;      % 残差
    df = 1 - e*cos(E);         % 导数
    dE = -f / df;              % 牛顿步长
    E = E + dE;
    
    if abs(dE) < tol || abs(f) < tol
        break;
    end
end
if iter == maxIter
    warning('开普勒方程迭代未收敛');
end

% 2. 计算轨道平面内的位置和速度 (近焦点坐标系)
% 位置分量 (x,y) 在轨道平面内，x指向近地点，y沿运动方向
r_orb = [a*(cos(E) - e);
         a*sqrt(1 - e^2)*sin(E);
         0];

r_norm = norm(r_orb);      % 或 a*(1 - e*cos(E))

% 速度分量
C = sqrt(mu*a) / r_norm;
v_orb = C * [-sin(E);
             sqrt(1 - e^2)*cos(E);
             0];

% 3. 构建从轨道平面到惯性系的旋转基向量 (P, Q, W)
% P: 指向近地点
% Q: 轨道平面内垂直于P，沿运动方向
% W: 角动量方向 (与轨道面法向)
P = [ cos(omega)*cos(RAAN) - sin(omega)*sin(RAAN)*cos(i);
      cos(omega)*sin(RAAN) + sin(omega)*cos(RAAN)*cos(i);
      sin(omega)*sin(i) ];

Q = [ -sin(omega)*cos(RAAN) - cos(omega)*sin(RAAN)*cos(i);
      -sin(omega)*sin(RAAN) + cos(omega)*cos(RAAN)*cos(i);
      cos(omega)*sin(i) ];

W = [ sin(RAAN)*sin(i);
     -cos(RAAN)*sin(i);
      cos(i) ];

% 4. 转换到惯性系
r = P*r_orb(1) + Q*r_orb(2);   % 第三分量为0
v = P*v_orb(1) + Q*v_orb(2);

end