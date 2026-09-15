function altitude = mission_minimum_altitude(solution)
% 端点、积分节点及所有定位到的径向极值；单位 m。
global R_E
states=solution.y(1:3,:);
if isfield(solution,'ye') && ~isempty(solution.ye), states=[states,solution.ye(1:3,:)]; end
altitude=min(vecnorm(states))-R_E;
end
