function fdt = Guoaccel(rho,u,du,cs2)
%Guo force term calculator
%Guo forcing 기반으로 추가 모멘텀을 계산합니다
%계산된 텀은 (1-1/2tau)항이 추가되지 않으므로 Relaxation(BGK, MRT) 단계에서
%계산이 될 수 있도록 해야 합니다

arguments
    rho(:,:)
    u (:,:,2)
    du (:,:,2)
    cs2 = 1/3
end

persistent xdir ydir weight
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
        weight=gpuArray(reshape([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],1,1,[]));
    end

cddu=du(:,:,1).*xdir+du(:,:,2).*ydir;
uddu=du(:,:,1).*u(:,:,1)+du(:,:,2).*u(:,:,2);
cdu=u(:,:,1).*xdir+u(:,:,2).*ydir;

fdt=weight.*(cddu-uddu+cdu.*cddu/cs2)/cs2.*rho;


end