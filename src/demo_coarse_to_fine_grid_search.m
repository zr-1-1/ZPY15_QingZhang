%% 无需 ATK 或碎片数据的搜索流程演示（示例评分不是清除能力指标）
addpath(fileparts(mfilename('fullpath')));
d2r = pi/180;
bounds = [7000e3,7450e3; .01,.04; 98*d2r,98.18*d2r; ...
    73*d2r,78*d2r; 345*d2r,445*d2r; 0,2*pi];
coarseStep = [150e3,.015,.09*d2r,2.5*d2r,50*d2r,90*d2r];
target = [7187500,.02125,98.0675*d2r,76.125*d2r,382.5*d2r,337.5*d2r];
scoreFcn = @(oe) demo_score(oe,target,coarseStep);
result = coarse_to_fine_grid_search(scoreFcn,bounds,coarseStep, ...
    struct('TopK',3,'Levels',2));
disp('最优根数 [m, 无量纲, rad, rad, rad, rad]：');
disp(result.BestOE);
disp(struct2table(result.History));

function score = demo_score(oe,target,scale)
delta = oe-target;
delta(4:6)=mod(delta(4:6)+pi,2*pi)-pi;
score = -sum((delta./scale).^2);
end
