
% % ===============
% %----------------
% % 볼츠만 실험기
% % 작성일:26-07-15
% % 작성자: 전도현
% % 설명: 격자-볼츠만 방법론을 적용한 대류 및 확산 시뮬레이션
% % 참조: 이 파일의 목표는 원기둥 벤치마크
% %----------------
% % ===============


%======= 시스템 세팅========
clear;
clear functions;

%========parameters (Physical)
% note: 무차원으로 입력 요망
%--------Lattice
Verlength=120; 
Horlength=360; 
lendel=1; 

Time=30000; 
timedel = 1; 
savestep = 20;


%--------Material Properties
rho_f=1;

H0=1;

tau_nu=0.501;

wall_rho=0.25;

%++++++ Initialization +++++++
%------ Space Genesis --------
verstep=Verlength/lendel;
horstep=Horlength/lendel;

frhofld=gpuArray(zeros(verstep+2,horstep+2));
u1=gpuArray(zeros(verstep+2,horstep+2,2));
du1=gpuArray(zeros(verstep+2,horstep+2,2));

xdir=reshape([0,1,0,-1,0,1,-1,-1,1],1,1,[]);
ydir=reshape([0,0,1,0,-1,1,1,-1,-1],1,1,[]);
weight=reshape([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36],1,1,[]);
%Just in Case
%dir: x> ; ^y
%     ^y
%    7 3 6
%    4 1 2 x>
%    8 5 9
%      
[X,Y]=meshgrid(1:horstep+2,1:verstep+2);


%------ Boundary (Wall) setting (or other obstacles later)
walls=gpuArray(uint8(zeros(verstep+2,horstep+2)));% 0==false==flow
%Note: matrix "walls" intentionally designated as uint8, for further wall
walls([1,end],:)=0;
walls=logical(((X-120).^2+(Y-60).^2).^0.5<15);

%------ Inital/boundary condition designation

%Initial rho & u distribution


%
%frhofld(logical(walls))=wall_rho;
% frhofld(:,:)=0.1300+0.0012*rand([verstep+2,horstep+2]);
frhofld(:,:)=rho_f;
% frhofld(:,:)=frhofld+(wall_rho-rho_g).*(2-tanh((Y-1)/2.2)-tanh((402-Y)/2.2));
u1(:)=0.01*rand(size(u1),'like',u1);
u1(:,1,1)=0.03+rand(size(u1(:,1,1)),"like",u1)*0.05;




%-----Initial Dist. Func. Calculation

feq=ustareq(frhofld,u1);
fdt=ustareq(frhofld,u1+du1)-feq;

f=feq;


%$$$$$$$$$$%%temporal savedata recall%%$$$$$$$$$$%
% load("temp.mat",'f','fdt','g','gdt')

% Initial Visualization
%------Image & Data Acquisition Initiation (tiledlayout 기반)
fig = figure('Color', 'w'); % 배경색 흰색으로 설정
tlo = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
hMaintitle=title(tlo, sprintf('%d step',0)); % 전체 제목

%--- 1번 타일: Velocity Norm (u1)
nexttile;
hImg(1) = imagesc(gather(vecnorm(u1, 2, 3))); 
colormap(gca, 'turbo'); colorbar;
clim([0, 0.08]);
title('Momentum Velocity Field (u1)');
axis equal tight;

%--- 2번 타일: Density Field (f)
nexttile;
hImg(2) = imagesc(gather(frhofld)); 
colormap(gca, 'turbo'); colorbar;
clim([0, 2]);
title('Liquid Density (frhofld)');
axis equal tight;

%--- 3번 타일: Velocity Norm (u2/du1)
nexttile;
hImg(3) = imagesc(gather(vecnorm(du1/2, 2, 3))); 
colormap(gca, 'turbo'); colorbar;
clim([0, 0.15]);
title('Interfacial Sustain (du1)');
axis equal tight;

%--- 4번 타일: Density Field (g)
nexttile;
hImg(4) = imagesc(gather(frhofld)); 
colormap(gca, 'turbo'); colorbar;
clim([0, 1.0]);
title('Gas Density (grhofld)');
axis equal tight;

% 공통 라벨 설정 (tiledlayout의 장점)
xlabel(tlo, 'X Grid Index');
ylabel(tlo, 'Y Grid Index');

drawnow;

tic



for i = 1:Time


    %------Collision (Relaxation
    f=DistributionMRT(f,feq,tau_nu)+fdt;
    f(:,2,:)=feq(:,1,:)+f(:,2,:)-feq(:,2,:);
    f(:,end,:)=f(:,end-1,:);

    % streaming and BBwalls
    f=stream(f);

    f=BBwalls(f,walls);


    %----Macroscopic Update
    [frhofld, u1]=rawdist(f);
    %u1(:)=u1+du1/2;

    %---boundary application
    frhofld(logical(walls))=wall_rho;
    frhofld(:,1)=1;
    frhofld(:,end)=1;
    u1(:,1,1)=0.05;
    u1(:,1,2)=0;
    u1(:,end,:)=0;

    %------Equilibrium and force term calculation
    feq=ustareq(frhofld,u1);
    fdt(:)=(ustareq(frhofld,u1+du1)-feq);
    
    %------Imaging------
    if rem(i,savestep)==0
        set(hImg(1), 'CData', gather(vecnorm(u1+du1/2,2,3)));
        set(hImg(2), 'CData', gather(frhofld));
        set(hImg(3), 'CData', gather(vecnorm(du1,2,3)));
        set(hImg(4), 'CData', gather(frhofld));
        set(hMaintitle, 'String', sprintf('%d step',i))
        drawnow limitrate;
    end
end

sum(u1(1:(verstep/2+1),:,:),"all")/sum(u1((verstep/2+2):end,:,:),"all")