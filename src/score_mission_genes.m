function [score,mask,shipMasks] = score_mission_genes(genes,catalog)
% 21 整数基因：先按列存放 3x4 轨道编号，再按列存放 3x3 转移时隙。
orbit=reshape(genes(1:12),3,4); slot=reshape(genes(13:21),3,3);
shipMasks=false(catalog.DebrisCount,3);
for window=1:4
    incoming=ones(3,1); outgoing=ones(3,1);
    if window>1, incoming=slot(:,window-1); end
    if window<4, outgoing=slot(:,window); end
    indices=sub2ind(catalog.MaskSizes{window},orbit(:,window),incoming,outgoing);
    shipMasks=shipMasks | catalog.Masks{window}(:,indices);
end
mask=any(shipMasks,2); score=nnz(mask);
end
