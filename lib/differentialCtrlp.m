function [new_ctrlp, new_cfg] = differentialCtrlp(order, ctrlp, cfg)
if order == 0
    new_ctrlp = ctrlp;
    new_cfg = cfg;
else
    dim = size(ctrlp);
    new_ctrlp = zeros(dim(1)-1, dim(2));
    for i=1:dim(1)-1
        new_ctrlp(i,:) = (cfg.k-1) * (ctrlp(i+1,:) - ctrlp(i, :)) / (cfg.knot(i+cfg.k) - cfg.knot(i+1));
    end
    new_cfg.k = cfg.k - 1;
    new_cfg.m = cfg.m;
    new_cfg.nc = dim(1)-1;
    new_cfg.knot = cfg.knot(2:length(cfg.knot)-1);
    [new_ctrlp, new_cfg] = differentialCtrlp(order-1, new_ctrlp, new_cfg);
end