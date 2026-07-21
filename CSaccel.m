function du = CSaccel(rhofield,T,cs2)
%   Carnahan-Starling pseudopotential model based Force calculator/phase
%   Carnahan-Starling pseudopotential model을 바탕으로 셀에 걸리는 의사퍼텐셜 가속(입자간 인척력)을 구합니다
%   벽 접착력은 디버그를 위해 따로 구하지 않습니다
arguments
    rhofield (:,:) %밀도장
    T           %온도(Tc)로
    cs2 = 1/3
end
persistent xdir ydir weight g
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
        weight=gpuArray(reshape([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],1,1,[]));
        g=-1;
    end

    pstar=rhofield.*T.*(1+rhofield+rhofield.^2-rhofield.^3)./(1-rhofield).^3-rhofield.^2*0.25-rhofield*cs2;


    psi=(max(2*pstar/g/cs2,0).^0.5);
    % wallID=logical(walls);
    Fx = zeros(size(rhofield),'like',rhofield);  
    Fy = zeros(size(rhofield),'like',rhofield);
    %Gfield=T.*~wallID;
    
    for i = 2:9
        nbr_p = circshift(psi, [-ydir(i), xdir(i)]);
        %nbr_g = circshift(Gfield, [-ydir(i), xdir(i)]);
        Fx = Fx - g*psi.* weight(i) .* nbr_p .* xdir(i);
        Fy = Fy - g*psi.* weight(i) .* nbr_p .* ydir(i);
    end
    
    du=zeros([size(rhofield),2],"like",rhofield);
    du(:,:,1)=-Fx./rhofield;
    du(:,:,2)=-Fy./rhofield;
    
end