function transfer = solve_j2_transfer(startState,targetState,duration)
% 两个零圈 Lambert 分支作初猜，在 model_rv 下微分修正到目标位置。
% 位置不可跳变；第二脉冲只匹配目标速度。全程高度至少 200 km。
global mu
validateattributes(startState,{'double'},{'real','finite','numel',6});
validateattributes(targetState,{'double'},{'real','finite','numel',6});
validateattributes(duration,{'double'},{'scalar','positive','finite'});
startState=startState(:); targetState=targetState(:);
transfer=struct('Feasible',false,'DeltaV',nan(3,2),'TotalDeltaV',inf, ...
    'PositionResidual',inf,'MinAltitude',-inf,'Branch',0,'Solution',[]);
options=mission_ode_options();
for branch=[1,-1]
    [velocity,ok]=lambert_guess(startState(1:3),targetState(1:3),duration,mu,branch);
    if ~ok, continue; end
    for iteration=1:12
        sol=mission_propagate([0,duration],[startState(1:3);velocity],options);
        residual=sol.y(1:3,end)-targetState(1:3);
        if norm(residual)<0.01, break; end
        jacobian=zeros(3);
        for axis=1:3
            perturbed=velocity; perturbed(axis)=perturbed(axis)+0.05;
            trial=mission_propagate([0,duration],[startState(1:3);perturbed],options);
            jacobian(:,axis)=(trial.y(1:3,end)-sol.y(1:3,end))/0.05;
        end
        if rcond(jacobian)<1e-10, break; end
        correction=jacobian\residual;
        if norm(correction)>3000, break; end
        velocity=velocity-correction;
    end
    % Reintegrate after the last correction, including iteration-limit exits.
    sol=mission_propagate([0,duration],[startState(1:3);velocity],options);
    error=norm(sol.y(1:3,end)-targetState(1:3));
    altitude=mission_minimum_altitude(sol);
    dv=[velocity-startState(4:6),targetState(4:6)-sol.y(4:6,end)];
    cost=sum(vecnorm(dv));
    if error<0.01 && altitude>=200e3 && cost<transfer.TotalDeltaV
        transfer=struct('Feasible',true,'DeltaV',dv,'TotalDeltaV',cost, ...
            'PositionResidual',error,'MinAltitude',altitude,'Branch',branch,'Solution',sol);
    end
end
end

function [velocity,ok] = lambert_guess(r1,r2,dt,mu,branch)
velocity=nan(3,1); ok=false;
n1=norm(r1); n2=norm(r2); cosine=max(-1,min(1,dot(r1,r2)/(n1*n2)));
A=branch*sqrt(n1*n2*(1+cosine));
if abs(A)<1e-8 || 1-cosine<1e-12, return; end
% Universal-variable zero-revolution time of flight, both transfer angles.
zs=linspace(-4*pi^2,4*pi^2-1e-5,600);
previous=nan; previousZ=nan;
for z=zs
    [time,~]=flight_time(z,A,n1,n2,mu);
    value=time-dt;
    if isfinite(value)
        if isfinite(previous) && previous*value<=0
            root=fzero(@(x) flight_time(x,A,n1,n2,mu)-dt,[previousZ,z]);
            [~,y]=flight_time(root,A,n1,n2,mu);
            f=1-y/n1; g=A*sqrt(y/mu);
            velocity=(r2-f*r1)/g; ok=all(isfinite(velocity)); return;
        end
        previous=value; previousZ=z;
    end
end
end

function [time,y] = flight_time(z,A,n1,n2,mu)
if abs(z)<1e-7
    C=1/2-z/24+z^2/720; S=1/6-z/120+z^2/5040;
elseif z>0
    q=sqrt(z); C=(1-cos(q))/z; S=(q-sin(q))/q^3;
else
    q=sqrt(-z); C=(cosh(q)-1)/(-z); S=(sinh(q)-q)/q^3;
end
y=n1+n2+A*(z*S-1)/sqrt(C);
if y<0 || C<=0, time=nan; return; end
time=((y/C)^1.5*S+A*sqrt(y))/sqrt(mu);
end
