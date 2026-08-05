function [grad_x,grad_y]=stencil_grad(rho)
arguments
    rho (:,:)
end
persistent xdir ydir weight
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
        weight=gpuArray(reshape([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],1,1,[]));
    end

    grad_x = zeros(size(rho),'like',rho);  
    grad_y = zeros(size(rho),'like',rho);
    %Gfield=T.*~wallID;
    
    for i = 2:9
        nbr_p = circshift(rho, [-ydir(i), xdir(i)]);
        grad_x = grad_x - weight(i) .* nbr_p .* xdir(i);
        grad_y = grad_y - weight(i) .* nbr_p .* ydir(i);
    end
    
    grad_x=grad_x.*3;
    grad_y=grad_y.*3;
    
end