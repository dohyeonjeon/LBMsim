function resr = stream(Matx)
%   Streaming step calculator
arguments
    Matx (:,:,9)
end
persistent xdir ydir fnew
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
        fnew=zeros(size(Matx),'like',Matx);
    end
    %dir: x> ; ^y
    %     ^y
    %    7 3 6
    %    4 1 2 x>
    %    8 5 9
    %      

    %WIP, change if boundary conditions are more complicated
    
    for i = 1:9
        fnew(:,:,i)=circshift(Matx(:,:,i),[-ydir(i),xdir(i)]);
    end
   
    resr=fnew;
end