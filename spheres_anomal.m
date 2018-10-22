clear
clc
h = 90;
prefix = {'/home/pzhang/chen/permeable_bed/'}
name = strcat(prefix,'spheres',num2str(h),'.txt');
spheres = load(char(name));
pos = spheres(2:end,:);
R = spheres(1,2);
N = spheres(1,3);
points = [3,59,74];
for i = 1:N
    dist(i,1) = sqrt(sum((pos(i,:) - points).^2));
    
end

[dd,I] = sort(dist);
for i = 1:N
    ii = I(i);
    ellipsoid(pos(ii,1),pos(ii,2),pos(ii,3),R,R,R);
    hold on
    axis equal
    drawnow
end
% shading interp
re = [0 1 0];
colormap(re)
alpha(0.5)
plot3(points(1),points(2),points(3),'r*')
