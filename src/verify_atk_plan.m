function verification=verify_atk_plan(plan,debris,c,options,solver)
% Integrate in ATK Inertial; compare with the existing J2000 debris library.
if nargin<5,solver=@ode113;end
inertialPlan=plan;
for ship=1:numel(plan.Ships)
    R=rotation(c,0); x=plan.Ships(ship).InitialState;
    inertialPlan.Ships(ship).InitialState=[R*x(1:3);R*x(4:6)];
    for k=1:size(plan.Ships(ship).Pulses,1)
        t=plan.Ships(ship).Pulses(k,1); R=rotation(c,t);
        inertialPlan.Ships(ship).Pulses(k,2:4)=(R*plan.Ships(ship).Pulses(k,2:4)')';
    end
end
verification=verify_mission_plan(inertialPlan,debris,options, ...
    @(t,y,o)solver(@(t,y)atk_calibrated_dynamics(t,y,c),t,y,o), ...
    @(t,y)to_j2000(t,y,c));
verification.Method='Empirical ATK inertial J2 replay, J2000 debris comparison; native ATK not verified';
end
function R=rotation(c,t)
w=c.FramePolynomial*[1;t/86400]; S=[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
R=eye(3)+S+.5*S*S;
end
function y=to_j2000(t,x,c)
y=x;
for k=1:numel(t)
    R=rotation(c,t(k)); y(:,k)=[R'*x(1:3,k);R'*x(4:6,k)];
end
end
