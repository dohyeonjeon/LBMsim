function del2=stencil_del2(rho)
arguments
    rho (:,:)
end
persistent xdir ydir weight
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
        weight=gpuArray(reshape([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],1,1,[]));
    end

    del2 = zeros(size(rho),'like',rho);

    for i = 1:numel(weight)
        del2 = del2 + weight(i) * (circshift(rho, [-ydir(i), xdir(i)]) - rho);
    end

    del2=del2*3/2;