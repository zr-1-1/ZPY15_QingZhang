function test_si_units
% 小规模 SI/弧度回归检查，不连接 ATK，不写入数据文件。
    global mu J_2 R_E
    oldValues = {mu, J_2, R_E};
    oldPath = path;
    cleanup = onCleanup(@() restoreState(oldValues, oldPath)); %#ok<NASGU>
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(projectRoot, 'src'));
    mu = 3.986004418e14;
    J_2 = -sqrt(5)*(-4.841653717360e-4);
    R_E = 6378137;

    % 圆轨道的解析状态检查米制尺度与 pi/2 弧度方向。
    a = 7000000;
    speed = sqrt(mu/a);
    state = orb_elements2rv([a, 0, 0, 0, 0, pi/2]);
    assert(norm(state(1:3)-[0,a,0]) < 1e-8);
    assert(norm(state(4:6)-[-speed,0,0]) < 1e-8);

    % 对实际 CSV 状态，用能量和角动量独立核对 a、e、i 的尺度。
    T = readtable(fullfile(projectRoot, 'data', 'debris_orbits.csv'));
    for k = 1:height(T)
        r = T{k,7:9}; v = T{k,10:12};
        oe = rv2coe(r,v);
        energy = dot(v,v)/2-mu/norm(r);
        h = cross(r,v);
        expectedA = -mu/(2*energy);
        expectedE = sqrt(1+2*energy*dot(h,h)/mu^2);
        assert(abs(oe(1)-expectedA) < 1e-7);
        assert(abs(oe(2)-expectedE) < 1e-10);
        assert(abs(cos(oe(3))-h(3)/norm(h)) < 1e-12);
        assert(all(isfinite(oe)) && all(oe(3:6)>=0) && all(oe(3:6)<2*pi));
    end

    % 无摄动时，传播一个开普勒周期恰好增加 2*pi 弧度。
    E0 = [a,.01,.5,.3,.2,1];
    savedJ2 = J_2;
    J_2 = 0;
    period = 2*pi*sqrt(a^3/mu);
    result = OE_scl_ptb(E0,period);
    assert(norm(result-[E0(1:5),E0(6)+2*pi]) < 1e-12);
    J_2 = savedJ2;
    assert(isequal(OE_scl_ptb(E0,0),E0));

    % 采用 p 与周期形式独立核对一阶 J2 平均角速率。
    p = a*(1-E0(2)^2);
    rate = (3*pi*J_2/period)*(R_E/p)^2;
    expectedRates = [-rate*cos(E0(3)), ...
        rate/2*(5*cos(E0(3))^2-1), ...
        2*pi/period+rate/2*sqrt(1-E0(2)^2)*(3*cos(E0(3))^2-1)];
    result = OE_scl_ptb(E0,86400);
    assert(max(abs(result(4:6)-(E0(4:6)+expectedRates*86400))) < 1e-12);

    % 所有长度放大 s、mu 放大 s^3，物理时间与角变化应保持一致。
    reference = result;
    scale = 2;
    mu = mu*scale^3;
    R_E = R_E*scale;
    scaledE0 = E0; scaledE0(1) = E0(1)*scale;
    result = OE_scl_ptb(scaledE0,86400);
    assert(max(abs(result(2:6)-reference(2:6))) < 1e-12);
    assert(result(1)==reference(1)*scale);

    % mu 缩小四倍且传播时间加倍，角变化应保持一致。
    % 这可检出摄动速率被额外乘一次平均运动的量纲错误。
    mu = mu/scale^3/4;
    R_E = R_E/scale;
    result = OE_scl_ptb(E0,2*86400);
    assert(max(abs(result-reference)) < 1e-12);
    fprintf('SI/radian checks passed (%d CSV states).\n',height(T));
end

function restoreState(oldValues, oldPath)
    global mu J_2 R_E
    mu = oldValues{1}; J_2 = oldValues{2}; R_E = oldValues{3};
    path(oldPath);
end
