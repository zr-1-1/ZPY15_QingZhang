function diagnose_reference_seeds()
% 仅核对论文表 1 的公开根数与本数据，不自动混入生产库。
% 表 1 最后一列角为 M（平近点角），已核对 PDF 第 6 页图像。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
% [论文航天器, 窗口, 原起点相对任务起点/h, a/km,e,i/deg,Omega/deg,omega/deg,M/deg]
raw=[1,1,0,7144,.009,98.10,75.80,40,63.7; ...
1,2,6,7178,.021,98.05,75.69,40,42.9; ...
1,3,12,7136,.016,98.05,76.03,20,34.1; ...
1,4,16,7242,.022,98.05,76.19,40,281.4; ...
2,1,0,7258,.030,98.05,75.45,50,136.1; ...
2,2,6,7260,.028,98.05,75.79,40,314.3; ...
2,3,9,7226,.026,98.05,75.91,40,52.1; ...
2,4,16,7266,.030,98.05,76.19,60,185.0; ...
3,1,0,7144,.015,98.20,75.80,20,285.8; ...
3,2,4,7214,.024,98.05,75.71,40,358.0; ...
3,3,12,7226,.022,98.05,76.03,60,330.7; ...
3,4,16,7232,.020,98.05,76.19,60,128.1];
windows=[0,8;6,14;12,20;16,24]*3600;
loaded=load(fullfile(root,'data','x_rv.mat'),'x_rv');
seeds=struct('PaperMother',{},'Window',{},'OriginalOE',{},'OriginalEpochSeconds',{}, ...
    'OE',{},'Score',{},'Intervals',{},'Interpretation',{});
for interpretation=1:2
for k=1:size(raw,1)
    window=raw(k,2); oe=raw(k,4:9).*[1000,1,pi/180,pi/180,pi/180,pi/180];
    % 正文称网格参数为真近点角，表头写 M；两种解释均作诊断，不能静默混用。
    if interpretation==2
        eccentric=2*atan2(sqrt(1-oe(2))*sin(oe(6)/2),sqrt(1+oe(2))*cos(oe(6)/2));
        oe(6)=mod(eccentric-oe(2)*sin(eccentric),2*pi);
    end
    original=oe; offset=windows(window,1)-raw(k,3)*3600;
    if offset>0
        sol=ode45(@model_rv,[0,offset],orb_elements2rv(oe)',mission_ode_options());
        oe=rv2coe(sol.y(1:3,end),sol.y(4:6,end));
    end
    sol=ode45(@model_rv,[0,diff(windows(window,:))],orb_elements2rv(oe)',mission_ode_options());
    intervals=mission_hit_intervals(sol,loaded.x_rv(:,:,windows(window,1)+1:windows(window,2)+1),windows(window,1));
    score=numel(unique(intervals(:,1)));
    seeds(end+1)=struct('PaperMother',raw(k,1),'Window',window,'OriginalOE',original, ...
        'OriginalEpochSeconds',raw(k,3)*3600,'OE',oe,'Score',score,'Intervals',intervals, ...
        'Interpretation',interpretation); %#ok<AGROW>
    fprintf('角度解释 %d，论文母船 %d 窗口 %d，当前标准窗口清除数 %d。\n',interpretation,raw(k,1),window,score);
end
end
save(fullfile(root,'data','mission','reference_seed_diagnostic.mat'),'seeds','raw','windows');
end
