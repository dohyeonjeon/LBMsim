function feq = ustareq(rho,ustar,cs2)
%충격량-전치된 유속 기반 균형 분포 설정
%   Guo forcing term에서는 기존 단계의 속도(f와 ei의 곱 총합 나누기 밀도) 대신
%   충격량을 받아 변동된 속도, u star를 사용합니다. 이를 기반으로 균형함수를 다시 계산하는 함수입니다.
arguments
    rho (:,:) %밀도장
    ustar (:,:,2) %속도장
    cs2=1/3 %격자음속^2
end
persistent xdir ydir weight
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
        weight=gpuArray(reshape([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],1,1,[]));
    end
    
    cdu=ustar(:,:,1).*xdir+ustar(:,:,2).*ydir;
    udu=ustar(:,:,1).^2+ustar(:,:,2).^2;

    feq=rho.*weight.*(1+cdu./cs2+cdu.^2./cs2.^2./2-udu./cs2./2);
end