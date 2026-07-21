function fout = DistTRT(f, feq, tau)
%MRT가 안 풀릴 때 TRT로 해 봅시다
%   MRT와 달리 두 개의 변수를 이용하는 단순화된 이완방정식입니다. 효율적이죠.
arguments
    f (:,:,9)
    feq (:,:,9)
    tau =1
end
persistent opp tao
    if isempty(opp)
    opp=gpuArray([1,4,5,2,3,8,9,6,7]);
    tao=4/tau;
    end
fopp=f(:,:,opp);
eopp=feq(:,:,opp);
fout = f +(feq-f)/tau+(eopp - fopp)/tao;
end