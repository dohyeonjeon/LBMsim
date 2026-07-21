function fnew = BBwalls(Matx, ID)
%   Bounce-back boundary applier
arguments
    Matx (:,:,9)
    ID (:,:) uint8
end
persistent xdir ydir opp
    if isempty(xdir)
        xdir=gpuArray([0,1,0,-1,0,1,-1,-1,1]);
        ydir=gpuArray([0,0,1,0,-1,1,1,-1,-1]);
        opp=gpuArray([1,4,5,2,3,8,9,6,7]);
    end

    %dir: x> ; ^y
    %     ^y
    %    7 3 6
    %    4 1 2 x>
    %    8 5 9
    %      

    
    fnew=zeros(size(Matx),'like',Matx);
    %tempf=zeros(size(Matx));
    wcanv=logical(ID);
    % fnew=tempf;
    for i = 1:9
        fnew(:,:,i)=Matx(:,:,i).*circshift(~wcanv,[-ydir(i),xdir(i)])+circshift(Matx(:,:,opp(i)).*wcanv,[-ydir(i),xdir(i)]);
    end

end