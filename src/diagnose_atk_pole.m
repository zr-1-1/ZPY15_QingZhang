function diagnose_atk_pole()
% Diagnostic only: estimate a fixed gravity pole from saved ATK coast endpoints.
load('data/atk_correction/dynamics_audit.mat','raw','rotation');
arcs=raw(cellfun(@(x)strcmp(x.kind,'CMCSPropagate'),raw));
first=arcs([1,8,15]);
options=optimset('Display','off','TolX',1e-11,'TolFun',1e-9,'MaxIter',150);
pole=fminsearch(@(p) objective(p,first,rotation),[0.003;0],options);
fprintf('Fitted inertial gravity pole x,y: %.15g %.15g\n',pole);
errors=zeros(numel(arcs),2);
for k=1:numel(arcs)
    e=arc_error(pole,arcs{k},rotation);
    errors(k,:)=[norm(e(1:3)),norm(e(4:6))];
end
disp(errors);
save('data/atk_correction/pole_diagnostic.mat','pole','errors');
end
function f=objective(p,arcs,R)
f=0;
for k=1:numel(arcs), e=arc_error(p,arcs{k},R); f=f+sum(e(1:3).^2); end
end
function e=arc_error(p,a,R)
x=a.InitialState.rv;
if strcmp(a.InitialState.frame,'J2000'),x=[R*x(1:3);R*x(4:6)];end
t0=datetime(a.InitialState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
t1=datetime(a.FinalState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
sol=ode113(@(t,y) dynamics(y,p),[0,seconds(t1-t0)],x, ...
    odeset('RelTol',2e-12,'AbsTol',1e-8,'MaxStep',60));
e=sol.y(:,end)-a.FinalState.rv;
end
function dy=dynamics(y,p)
r=y(1:3); v=y(4:6); n=norm(r); pole=[p(:);sqrt(1-sum(p.^2))];
z=dot(r,pole); mu=3.986004415e14; radius=6378136.3; j2=sqrt(5)*4.841653717360e-4;
acc=-mu/n^3*r+1.5*mu*j2*radius^2/n^5*((5*z^2/n^2-1)*r-2*z*pole);
dy=[v;acc];
end
