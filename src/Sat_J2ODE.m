function dState = Sat_J2ODE(t, State)
% J2摄动下航天器受力的微分方程
    C20 = -4.841653717360e-04;
    J_2 = -sqrt(2*2+1)*C20; % 由完全归一化球谐系数得到J2摄动系数，此处为正数
    R_E_m = 6378137; % 地球赤道平均半径，单位：m
    mu = 3.986004418e14; % 地球引力常数，单位：(m^3)/(s^2)

    r = norm(State(1:3));
    pos = State(1:3);
    vel = State(4:6);
    a_central = - mu/r^3.*pos;
    C = 1.5*mu*J_2*R_E_m^2/r^5;
    x = pos(1); y = pos(2); z = pos(3);
    a_J2 = C.*[x*(5*(z/r)^2-1); y*(5*(z/r)^2-1); z*(5*(z/r)^2-3)];
    dState = [vel; a_central+a_J2];

end