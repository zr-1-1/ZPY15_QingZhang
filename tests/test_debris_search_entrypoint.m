function test_debris_search_entrypoint
% 在临时目录运行独立脚本，验证读取、时间截取、保存与元数据；不覆盖真实结果。
repository=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repository,'src'));
fixtureRoot=tempname;
mkdir(fixtureRoot); mkdir(fullfile(fixtureRoot,'src')); mkdir(fullfile(fixtureRoot,'data'));
cleanup=onCleanup(@()remove_fixture(fixtureRoot)); %#ok<NASGU>
global J_2 R_E mu
J_2=-sqrt(5)*(-4.841653717360e-4); R_E=6378137; mu=3.986004418e14;
oe=[7200e3,.025,1.71,1.3,.4,.6];
mother=orb_elements2rv(oe)';
% 数据第 3 帧对应母航天器初始时刻，保证两个历元没有混淆。
x_rv=zeros(6,8,5);
sol=ode45(@model_rv,[0,2],mother, ...
    odeset('RelTol',1e-9,'AbsTol',[1e-3;1e-3;1e-3;1e-6;1e-6;1e-6]));
for frame=3:5
    x_rv(:,:,frame)=repmat(deval(sol,frame-3),1,8);
end
save(fullfile(fixtureRoot,'data','x_rv.mat'),'x_rv');
source=fileread(fullfile(repository,'src','run_debris_grid_search.m'));
injection=sprintf(['searchStartUTC=stateDataStartUTC+seconds(2);\n' ...
    'searchEndUTC=stateDataStartUTC+seconds(4);\n' ...
    'searchBounds=[7200e3,.025,1.71,1.3,.4,.6]''*[1,1];\n' ...
    'executionOptions.EnableParallel=false;\n' ...
    'searchOptions.Verbose=false;\n']);
source=strrep(source,'stateFile = fullfile',[injection,'stateFile = fullfile']);
script=fullfile(fixtureRoot,'src','run_debris_grid_search.m');
fid=fopen(script,'w','n','UTF-8'); assert(fid>=0);
fprintf(fid,'%s',source); fclose(fid);
run(script);
saved=load(fullfile(fixtureRoot,'data','grid_search_result.mat'),'searchResult');
library=load(fullfile(fixtureRoot,'data','preferred_orbit_library.mat'),'preferredOrbitLibrary');
r=saved.searchResult;
assert(r.BestScore==8 && numel(r.History)==3);
assert(r.Evaluations==1 && r.CacheHits==2);
assert(r.OrbitEpochUTC==r.StateDataStartUTC+seconds(2));
assert(r.SearchEndUTC==r.StateDataStartUTC+seconds(4));
assert(all(r.FirstClearTime==0) && all(r.FirstClearUTC==r.OrbitEpochUTC));
assert(isequal(library.preferredOrbitLibrary.Orbits,r.PreferredOrbits));
assert(library.preferredOrbitLibrary.SearchStartUTC==r.SearchStartUTC);
assert(library.preferredOrbitLibrary.SearchEndUTC==r.SearchEndUTC);
assert(~r.Execution.ParallelUsed);
fprintf('test_debris_search_entrypoint passed\n');
end

function remove_fixture(folder)
% 仅删除本测试已知文件，然后移除空目录；不递归删除计算路径。
files={fullfile(folder,'src','run_debris_grid_search.m'), ...
    fullfile(folder,'data','x_rv.mat'), ...
    fullfile(folder,'data','grid_search_result.mat'), ...
    fullfile(folder,'data','preferred_orbit_library.mat')};
for k=1:numel(files)
    if isfile(files{k}), delete(files{k}); end
end
if contains([path,pathsep],[fullfile(folder,'src'),pathsep])
    rmpath(fullfile(folder,'src'));
end
rmdir(fullfile(folder,'src')); rmdir(fullfile(folder,'data')); rmdir(folder);
end
