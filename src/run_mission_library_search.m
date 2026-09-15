function run_mission_library_search(levels)
% 补齐 README 三个缺失窗口，不覆盖已有单窗口结果。支持按窗口恢复。
% levels=0 用于粗搜诊断；正式运行 levels=2，严格 >5 细化、>=8 入库。
if nargin<1, levels=2; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
output=fullfile(root,'data','mission'); if ~isfolder(output), mkdir(output); end
diary(fullfile(output,sprintf('library_search_L%d.log',levels)));
finishDiary=onCleanup(@()diary('off')); %#ok<NASGU>
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
physics=struct('J_2',J_2,'R_E',R_E,'mu',mu);
epoch=datetime(2030,11,14,8,0,0,'TimeZone','UTC');
windows=[0,8;6,14;12,20;16,24]*3600;
d=pi/180;
bounds=[7000e3,7450e3;.01,.04;98*d,98.18*d;73*d,78*d;345*d,445*d;0,2*pi];
step=[150e3,.015,.09*d,2.5*d,50*d,90*d];
opts=struct('TopK',3,'Levels',levels,'RefineFactor',2,'MaxEvaluations',200000, ...
    'Verbose',true,'RefineThreshold',5,'ArchiveThreshold',8,'BatchSize',128, ...
    'UseCache',true,'ProgressIntervalSeconds',10);
execution=struct('EnableParallel',true,'NumWorkers',6,'WorkerThreads',1, ...
    'ChunkSize',512,'ResizeExistingPool',false);
loaded=load(fullfile(root,'data','x_rv.mat'),'x_rv');
assert(isequal(size(loaded.x_rv),[6,345,86401]),'需要 345 个碎片、24 h、1 s 轨迹。');
csv=readtable(fullfile(root,'data','debris_orbits.csv'));
expected=[csv.PositionX,csv.PositionY,csv.PositionZ,csv.VelocityX,csv.VelocityY,csv.VelocityZ]';
assert(max(abs(loaded.x_rv(:,:,1)-expected),[],'all')<.1,'轨迹首帧与 CSV 不符，需确认数据来源。');
assert(isequal(csv.DebrisID(:),(1:345)'),'碎片编号顺序不一致。');
for window=2:4
    filename=fullfile(output,sprintf('window%d_L%d.mat',window,levels));
    if isfile(filename), fprintf('保留已完成窗口：%s\n',filename); continue; end
    fprintf('开始窗口 %d：%s 至 %s UTC\n',window,char(epoch+seconds(windows(window,1))),char(epoch+seconds(windows(window,2))));
    frames=windows(window,1)+1:windows(window,2)+1;
    coarseFile=fullfile(output,sprintf('window%d_L0.mat',window));
    if levels>0 && isfile(coarseFile)
        coarse=load(coarseFile,'result');
    else
        coarse=struct();
    end
    if isfield(coarse,'result') && coarse.result.BestScore<=5 && ...
            isequal(coarse.result.Bounds,bounds) && isequal(coarse.result.History(1).Step,step)
        result=coarse.result; result.Options.Levels=levels;
        result.StopReason='NoEligibleSeeds';
        fprintf('复用已完成粗搜；最高分 <=5，按门槛停止细搜。\n');
    else
        result=execute_debris_grid_search(loaded.x_rv(:,:,frames),(0:numel(frames)-1)', ...
            bounds,step,opts,execution,physics);
    end
    result.OrbitEpochUTC=epoch+seconds(windows(window,1));
    result.SearchStartUTC=result.OrbitEpochUTC; result.SearchEndUTC=epoch+seconds(windows(window,2));
    result.StateDataStartUTC=epoch; result.StateStepSeconds=1;
    result.TargetMin=5; result.TargetNum=8; result.Physics=physics;
    save(filename,'result');
    fprintf('窗口 %d 完成，最佳 %d，入库 %d 条。\n',window,result.BestScore,size(result.PreferredOrbits,1));
end
end
