function audit_atk_dynamics()
% Compare each saved ATK arc independently; never reset states in a mission replay.
addpath('src');
global mu R_E J_2
raw=jsondecode(fileread('data/atk_correction/original_segments.json'));
if isstruct(raw), raw=num2cell(raw); end
burns=raw(cellfun(@(x)strcmp(x.kind,'CMCSManeuver'),raw));
B=zeros(3,numel(burns)); A=B;
for k=1:numel(burns)
    b=burns{k}; B(:,k)=b.burn;
    A(:,k)=b.FinalState.rv(4:6)-b.InitialState.rv(4:6);
end
[u,~,v]=svd(A*B'); rotation=u*v';
fprintf('J2000 -> Inertial rotation:\n'); disp(rotation);
fprintf('All-burn rotation fit max residual %.9g m/s\n',max(vecnorm(A-rotation*B)));
arcs=raw(cellfun(@(x)strcmp(x.kind,'CMCSPropagate'),raw));
parameters=[3.986004418e14,6378137;3.986004415e14,6378136.3];
results=[];
for model=1:2
    mu=parameters(model,1); R_E=parameters(model,2); J_2=-sqrt(5)*(-4.841653717360e-4);
    for k=1:numel(arcs)
        a=arcs{k}; x=a.InitialState.rv;
        if strcmp(a.InitialState.frame,'J2000'), x=[rotation*x(1:3);rotation*x(4:6)]; end
        t0=datetime(a.InitialState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
        t1=datetime(a.FinalState.utc,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
        sol=mission_propagate([0,seconds(t1-t0)],x,mission_ode_options());
        error=sol.y(:,end)-a.FinalState.rv;
        results(end+1,:)=[model,a.mother,k,seconds(t1-t0),norm(error(1:3)),norm(error(4:6))]; %#ok<AGROW>
    end
end
resultTable=array2table(results,'VariableNames',{'Model','Mother','Arc','Duration_s','PositionError_m','VelocityError_mps'});
disp(resultTable);
writetable(resultTable,'data/atk_correction/arc_dynamics_audit.csv');
save('data/atk_correction/dynamics_audit.mat','rotation','resultTable','raw','parameters');
end
