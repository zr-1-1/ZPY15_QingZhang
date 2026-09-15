clc;
clear;
load("high_capture_segments.mat");
bestob = [7 5 1 3; 2 10 3 11; 8 2 5 14];
global J_2 R_E_m mu
J_2 = 1.082627e-3;  R_E_m = 6378137;  mu = 3.986004418e14;

%%
Lambert_cell = cell(3,4);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
window_st = [0, 6, 12, 16].*3600;
for sd = 1:3
    for td = 1:4
        oe = high_all{td}(bestob(sd,td)).oe;
        [r, v] = orb_elements2rv(oe,mu);
        ct = high_all{td}(bestob(sd,td)).capture_times;
        tspan = [0, [min(ct),max(ct)]-window_st(td)];
        [~, State] = ode78(@Sat_J2ODE, tspan, [r;v], opts);
        Sat_r = State(2:3,1:3)';
        s = struct('r', num2cell(Sat_r, 1), 't', {min(ct),max(ct)});
        Lambert_cell{sd,td} = s;
    end
end