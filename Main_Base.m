
% % ===============
% %----------------
% % 볼츠만 실험기
% % 작성일:26-07-14
% % 작성자: 전도현
% % 설명: 격자-볼츠만 방법론을 적용한 대류 및 확산 시뮬레이션
% % 참조: 이 폴더의 목표는 안정성 조건의 탐구와 확보임
% %----------------
% % ===============


%======= 시스템 세팅========
clear;
clear functions;

%========parameters (Physical)
% note: 무차원으로 입력 요망
%--------Lattice
Verlength=200; 
Horlength=200; 
lendel=1; 

Time=100; 
timedel = 1; 
savestep = 100;

Warmup_steps=50000;
warmup_Diff=1.5;

%--------Material Properties
rho_l=0.35955; rho_g=0.03015;
T=0.080*0.25;

H0=1;

tau_nu=0.6; tau_D=0.7;

wall_rho=0.25;
wall_rho_2=0.018;

rho_i=wall_rho*1.05
%++++++ Initialization +++++++
%------ Space Genesis --------
verstep=Verlength/lendel;
horstep=Horlength/lendel;

frhofld=gpuArray(zeros(verstep+2,horstep+2));
fphifld=gpuArray(zeros(verstep+2,horstep+2));
u1=gpuArray(zeros(verstep+2,horstep+2,2));
du1=gpuArray(zeros(verstep+2,horstep+2,2));

grhofld=gpuArray(zeros(verstep+2,horstep+2));
u2=gpuArray(zeros(verstep+2,horstep+2,2));
du2=gpuArray(zeros(verstep+2,horstep+2,2));

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
wall_mask=logical(((X-70).^2+(Y-90).^2).^0.5<20);
wall_mask(Y>107)=1;
walls(wall_mask)=1;

%------ Inital/boundary condition designation

%Initial rho & u distribution


%
%frhofld(logical(walls))=wall_rho;
% frhofld(:,:)=0.1300+0.0012*rand([verstep+2,horstep+2]);
frhofld(:,:)=(rho_g+(rho_l-rho_g)*(1-tanh((((X-100).^2+(Y-100).^2).^0.5-10)/4))/2);
frhofld(:,:)=(frhofld+(rho_l-rho_g)*(1-tanh((((X-70).^2+(Y-90).^2).^0.5-20)/4))/2);
% frhofld(:,:)=frhofld+(wall_rho-rho_g).*(2-tanh((Y-1)/2.2)-tanh((402-Y)/2.2));
u1(:)=CSaccel(frhofld,T);
du1(:)=0;


grhofld(:,:)=0.2;
grhofld(logical(walls))=wall_rho_2;
u2(:)=0;
du2(:)=0;



%-----Initial Dist. Func. Calculation

feq=ustareq(frhofld,u1);
fdt=ustareq(frhofld,u1+du1)-feq;

geq=ustareq(grhofld,u1);
gdt=ustareq(grhofld,u1+du1)-geq;

f=feq;
g=geq;


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
clim([0, 0.15]);
title('Momentum Velocity Field (u1)');
axis equal tight;

%--- 2번 타일: Density Field (f)
nexttile;
hImg(2) = imagesc(gather(frhofld)); 
colormap(gca, 'turbo'); colorbar;
clim([0, 0.5]);
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
hImg(4) = imagesc(gather(grhofld)); 
colormap(gca, 'turbo'); colorbar;
clim([0, 1.0]);
title('Gas Density (grhofld)');
axis equal tight;

% 공통 라벨 설정 (tiledlayout의 장점)
xlabel(tlo, 'X Grid Index');
ylabel(tlo, 'Y Grid Index');

drawnow;

tic

%+++++++ Warmup Loop +++++++

for i = 1:Warmup_steps


    %------Collision (Relaxation
    f=DistributionMRT(f,feq,warmup_Diff)+fdt;
    g=DistTRT(g,geq,warmup_Diff);

    % streaming and BBwalls
    f=stream(f);
    g=stream(g);

    f=BBwalls(f,walls);
    g=BBwalls(g,walls);


    %----Macroscopic Update
    [frhofld, u1]=rawdist(f);
    % [fphifld, ~]=rawdist(phi);
    du1=CSaccel(frhofld,T);
    %u1(:)=u1+du1/2;

    [grhofld, u2]=rawdist(g);
    du2(:)=0;
    %u2(:)=u2+du2/2;

    %---boundary application
    frhofld(logical(walls))=wall_rho;
    grhofld(frhofld>rho_i)=H0;
    grhofld([1,end],:)=0.2*H0;
    grhofld(:,[1,end])=0.2*H0;

    %------Equilibrium and force term calculation
    feq=ustareq(frhofld,u1);
    fdt(:)=(ustareq(frhofld,u1+du1)-feq);
    % phieq=ustareq(fphifld,u1);

    geq=ustareq(grhofld,u1+du1/2);
    gdt(:)=ustareq(grhofld,u1+du1)-geq;

    %------Imaging------
    if rem(i,savestep)==0
        set(hImg(1), 'CData', gather(vecnorm(u1+du1/2,2,3)));
        set(hImg(2), 'CData', gather(frhofld));
        set(hImg(3), 'CData', gather(vecnorm(du1,2,3)));
        set(hImg(4), 'CData', gather(grhofld));
        set(hMaintitle, 'String', sprintf('%d step',i))
        drawnow limitrate;
    end
end


for i = 1:Time


    %------Collision (Relaxation
    f=DistributionMRT(f,feq,tau_nu)+fdt;
    g=DistTRT(g,geq,tau_D);

    % streaming and BBwalls
    f=stream(f);
    g=stream(g);

    f=BBwalls(f,walls);
    g=BBwalls(g,walls);


    %----Macroscopic Update
    [frhofld, u1]=rawdist(f);
    % [fphifld, ~]=rawdist(phi);
    du1=CSaccel(frhofld,T);
    %u1(:)=u1+du1/2;

    [grhofld, u2]=rawdist(g);
    du2(:)=0;
    %u2(:)=u2+du2/2;

    %---boundary application
    frhofld(logical(walls))=wall_rho;
    grhofld(frhofld>rho_i)=H0;
    grhofld([1,end],:)=0.2*H0;
    grhofld(:,[1,end])=0.2*H0;

    %------Equilibrium and force term calculation
    feq=ustareq(frhofld,u1);
    fdt(:)=(ustareq(frhofld,u1+du1)-feq);
    % phieq=ustareq(fphifld,u1);

    geq=ustareq(grhofld,u1+du1/2);
    gdt(:)=ustareq(grhofld,u1+du1)-geq;

    %------Imaging------
    if rem(i,savestep)==0
        set(hImg(1), 'CData', gather(vecnorm(u1+du1/2,2,3)));
        set(hImg(2), 'CData', gather(frhofld));
        set(hImg(3), 'CData', gather(vecnorm(du1,2,3)));
        set(hImg(4), 'CData', gather(grhofld));
        set(hMaintitle, 'String', sprintf('%d step',i))
        drawnow limitrate;
    end
end