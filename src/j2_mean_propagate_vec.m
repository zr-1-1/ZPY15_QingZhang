function [r_eci, v_eci] = j2_mean_propagate_vec(oe0, t_array, mu, J_2, R_E_m)
% J2_MEAN_PROPAGATE_VEC  与 OE_scl_ptb 完全一致的 J2 平均根数递推（向量化）
% 输入:
%   oe0      : [a, e, i, Om, w, M]  (m, rad)
%   t_array  : 时间向量 (s)
%   mu, J_2, R_E_m : 引力常数、J2、赤道半径（与 Sat_J2ODE 保持一致）
% 输出:
%   r_eci, v_eci : 3×N (m, m/s)

a0 = oe0(1); e0 = oe0(2); i0 = oe0(3);
Om0 = oe0(4); w0 = oe0(5); M0 = oe0(6);

% --- 与 OE_scl_ptb 同一公式 ---
C_J2 = 1.5 * J_2 * R_E_m^2 * sqrt(mu) * a0^(-3.5);

Om_t = Om0 - C_J2 * cos(i0) / (1 - e0^2)^2 * t_array;
w_t  = w0  + C_J2 * (2 - 2.5*sin(i0)^2) / (1 - e0^2)^2 * t_array;
M_t  = M0  + sqrt(mu/a0^3) * ...
       (1 + 0.5*C_J2*(1 - e0^2)^(-1.5)*(3*cos(i0)^2 - 1)) * t_array;

% --- 解 Kepler 方程（向量化牛顿迭代） ---
E = M_t;
for k = 1:10
    E = E - (E - e0*sin(E) - M_t) ./ (1 - e0*cos(E));
end

cosE = cos(E); sinE = sin(E);
cosnu = (cosE - e0) ./ (1 - e0*cosE);
sinnu = sqrt(1 - e0^2) * sinE ./ (1 - e0*cosE);
r_mag = a0 * (1 - e0*cosE);

p0 = a0 * (1 - e0^2);
h  = sqrt(mu * p0);

% 近焦点系位置/速度
x_orb  = r_mag .* cosnu;
y_orb  = r_mag .* sinnu;
vx_orb = -h/p0 .* sinnu;
vy_orb =  h/p0 .* (e0 + cosnu);

% 旋转基（与 orb_elements2rv 中的 P、Q 完全一致）
cO = cos(Om_t); sO = sin(Om_t);
cw = cos(w_t);  sw = sin(w_t);
cI = cos(i0);   sI = sin(i0);

P11 =  cO.*cw - sO.*sw*cI;
P12 = -cO.*sw - sO.*cw*cI;
P21 =  sO.*cw + cO.*sw*cI;
P22 = -sO.*sw + cO.*cw*cI;
P31 =  sw*sI;
P32 =  cw*sI;

r_eci = [P11.*x_orb + P12.*y_orb;
         P21.*x_orb + P22.*y_orb;
         P31.*x_orb + P32.*y_orb];

v_eci = [P11.*vx_orb + P12.*vy_orb;
         P21.*vx_orb + P22.*vy_orb;
         P31.*vx_orb + P32.*vy_orb];
end