function distr=lapcalc(rho)
arguments
    rho (:,:)
end
weight=[16,4,4,4,4,1,1,1,1]/36;
xdir=[0,1,0,-1,0,1,-1,-1,1];
ydir=[0,0,1,0,-1,1,1,-1,-1];
distr=ustareq(rho,0);
target=zeros([size(rho)]);
for i = 1:9
    target(:,:) = weight(i) * circshift(distr(:,:,i),[-ydir(i),xdir(i)]);
end
distr=target;
end

