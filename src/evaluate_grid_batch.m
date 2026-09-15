function scores = evaluate_grid_batch(scoreFcn, candidates, useParallel, maxWorkers)
% 每行独立评分，结果严格按候选行号返回；异常直接中止，避免不完整结果。
scores=zeros(size(candidates,1),1);
if nargin<4, maxWorkers=4; end
if useParallel && ~isempty(candidates)
    parfor (row=1:size(candidates,1),maxWorkers)
        score=scoreFcn(candidates(row,:));
        validateattributes(score,{'numeric'},{'scalar','real','finite'});
        scores(row)=double(score);
    end
else
    for row=1:size(candidates,1)
        score=scoreFcn(candidates(row,:));
        validateattributes(score,{'numeric'},{'scalar','real','finite'});
        scores(row)=double(score);
    end
end
end
