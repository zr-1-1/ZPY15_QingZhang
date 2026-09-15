function solution = mission_propagate(times,state,options)
% 高阶 Adams 积分用于多次大脉冲后的开环回放，降低误差逐段放大。
% 动力学仍严格调用 model_rv；只更换数值积分器。
if nargin<3, options=mission_ode_options(); end
solution=ode113(@model_rv,times,state,options);
end
