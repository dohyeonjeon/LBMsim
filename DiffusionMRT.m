function resr = DiffusionMRT(disorig,diseq,tau)
%   MRT relaxation function for diffusion
%   3개 혹은 5개의 인수를 받습니다
%   1~3개의 인수는 각각 원래 분포, 균형 분포, tau입니다
%   힘 상수와 확산계수를 추가로 받으면 그 값도 계산에 넣습니다
%dir: x> ; ^y
%     ^y
%    7 3 6
%    4 1 2 x>
%    8 5 9
arguments
    disorig (:,:,9) 
    diseq (:,:,9)
    tau =1 ;
end

   
persistent M luM svec
    if isempty(M)
        M=gpuArray([ 1, 1, 1, 1, 1, 1, 1, 1, 1;    % m0 = density
           -4,-1,-1,-1,-1, 2, 2, 2, 2;    % m1 = e (energy)
            4,-2,-2,-2,-2, 1, 1, 1, 1;    % m2 = epsilon (energy²)
            0, 1, 0,-1, 0, 1,-1,-1, 1;    % m3 = jx (momentum x)
            0,-2, 0, 2, 0, 1,-1,-1, 1;    % m4 = qx (energy flux x)
            0, 0, 1, 0,-1, 1, 1,-1,-1;    % m5 = jy (momentum y)
            0, 0,-2, 0, 2, 1, 1,-1,-1;    % m6 = qy (energy flux y)
            0, 1,-1, 1,-1, 0, 0, 0, 0;    % m7 = pxx (normal stress)
            0, 0, 0, 0, 0, 1,-1, 1,-1]);   % m8 = pxy (shear stress)
        luM=M\eye(9);
        svec=gpuArray(diag([1.0, 1.0/tau, 1.0/tau, 1.0/tau, 1.0/tau, 1.0/tau, 1.0/tau, 1.0, 1.0]));
    end

svec([31,51])=1/tau;
morig=M*(reshape(disorig,[],9).');
mstar=morig-svec*(morig-M*reshape(diseq,[],9).');%+(eye(9)-svec/2)*(M*reshape(fdt,[],9).');
resr=reshape((luM*mstar).',size(disorig));

end