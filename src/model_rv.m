% 母航天器动力学模型
function [dxdt] = model_rv(t,x)
    % 导入全局变量，保持参数一致性
    global J_2 R_E w_E mu

    r0 = x(1:3,:);
    norm_r = norm(r0);
    v0 = x(4:6,:);

    g = -(mu./(norm_r)^3).*r0;
    a_J2 = [1 1 1]' .* (1.5*mu*J_2*(R_E^2)/(norm_r^5));
    a_J2(1,1) = a_J2(1,1).*(5*r0(1)*(r0(3)^2)/(norm_r^2)-r0(1));
    a_J2(2,1) = a_J2(2,1).*(5*r0(2)*(r0(3)^2)/(norm_r^2)-r0(2));
    a_J2(3,1) = a_J2(3,1).*(5*(r0(3)^3)/(norm_r^2)-3*r0(3));
    dxdt = [v0(1);v0(2);v0(3);g+a_J2];
end