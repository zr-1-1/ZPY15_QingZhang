function calibrate_atk_dynamics()
% Empirical diagnostic adapter, calibrated on Mother1/2 and held out on Mother3.
% This is not a substitute for a native ATK validation.
load('data/atk_correction/dynamics_audit.mat','raw');
load('data/atk_correction/pole_diagnostic.mat','pole');
epoch=datetime(2030,11,14,8,0,0);
burns=raw(cellfun(@(x)strcmp(x.kind,'CMCSManeuver'),raw));
design=[]; rhs=[];
for k=1:numel(burns)
    b=burns{k}; v=b.burn(:); t=seconds(datetime(b.InitialState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS')-epoch)/86400;
    skew=[0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
    design=[design;-skew,-t*skew]; %#ok<AGROW>
    rhs=[rhs;b.FinalState.rv(4:6)-b.InitialState.rv(4:6)-v]; %#ok<AGROW>
end
frame=reshape(design\rhs,3,2);
arcs=raw(cellfun(@(x)strcmp(x.kind,'CMCSPropagate'),raw));
for k=1:numel(arcs)
    a=arcs{k};
    a.t0=seconds(datetime(a.InitialState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS')-epoch);
    a.t1=seconds(datetime(a.FinalState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS')-epoch);
    a.x=a.InitialState.rv;
    if strcmp(a.InitialState.frame,'J2000')
        R=frame_rotation(frame,a.t0); a.x=[R*a.x(1:3);R*a.x(4:6)];
    end
    arcs{k}=a;
end
calibration=struct('PolePolynomial',[pole(:),zeros(2,2)],'FramePolynomial',frame, ...
    'Mu',3.986004415e14,'GravityRadius',6378136.3,'J2',sqrt(5)*4.841653717360e-4, ...
    'Status','Empirical ATK adapter; native rerun required');
train=arcs(cellfun(@(x)x.mother<3,arcs));
p=calibration.PolePolynomial(:);
for iteration=1:4
    e=residual(p,train,calibration); jac=zeros(numel(e),numel(p));
    for j=1:numel(p)
        q=p; q(j)=q(j)+1e-7; jac(:,j)=(residual(q,train,calibration)-e)/1e-7;
    end
    step=jac\e; p=p-step;
    fprintf('Calibration iteration %d, position RMS %.9g m\n',iteration,sqrt(mean(e.^2)));
end
calibration.PolePolynomial=reshape(p,2,3);
e=reshape(residual(p,arcs,calibration),3,[])';
calibration.ArcPositionErrors=vecnorm(e,2,2);
disp(calibration.PolePolynomial); disp(calibration.ArcPositionErrors');
save('data/atk_correction/calibration.mat','calibration','arcs');
end
function e=residual(p,arcs,c)
c.PolePolynomial=reshape(p,2,3); e=[];
for k=1:numel(arcs)
    a=arcs{k}; sol=ode113(@(t,y) atk_calibrated_dynamics(t,y,c),[a.t0,a.t1],a.x, ...
        odeset('RelTol',3e-14,'AbsTol',[1e-9;1e-9;1e-9;1e-12;1e-12;1e-12],'MaxStep',15));
    e=[e;sol.y(1:3,end)-a.FinalState.rv(1:3)]; %#ok<AGROW>
end
end
function R=frame_rotation(f,t)
w=f*[1;t/86400]; S=[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0]; R=expm(S);
end
