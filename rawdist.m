function [Rho, u] = rawdist(Matx)
%   Equilibrium distribution (Boltzmann dist.) finder
%%%%%%

arguments
    Matx (:,:,9)
end
persistent xdir ydir
    if isempty(xdir)
        xdir=gpuArray(reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]));
        ydir=gpuArray(reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]));
    end


    %dir: x> ; ^y
    %     ^y
    %    7 3 6
    %    4 1 2 x>
    %    8 5 9
    %      

    Rho=sum(Matx,3);
    Rho=max(Rho,eps);
    u_x=sum(Matx.*xdir,3)./Rho;
    v_y=sum(Matx.*ydir,3)./Rho;
    u=zeros([size(Rho),2],'like',Matx);
    u(:,:,1)=u_x;
    u(:,:,2)=v_y;

end
