function options = mission_ode_options()
% 高精度传播，同时定位径向速度零点以检查弧段最低高度。
options=odeset('RelTol',2e-13,'AbsTol',[1e-8;1e-8;1e-8;1e-11;1e-11;1e-11], ...
    'MaxStep',30,'Events',@radial_extremum);
end

function [value,terminal,direction] = radial_extremum(~,state)
value=dot(state(1:3),state(4:6));
terminal=0;
direction=0;
end
