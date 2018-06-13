clear
clc
prefix_name = {'/home/pzhang/chen/permeable_bed/'};
middle_name = {'test_carve_'};
posx = 200.5;
posy = 299.5;

num = 1:554;
name = strcat(prefix_name,middle_name,num2str(num(1),'%04d'),'.h5');
[data,domain] = getData(char(name),char('/Flbm_0'),2,true);
Nx = domain.Nx;
Ny = domain.Ny;
R = 20;
[x,y] = meshgrid(0:Nx-1, 0:Ny-1);
xx = x-posx;
yy = y-posy;
kkk = find((xx.^2+yy.^2)>(1.2*R^2));
%plot(x(kkk),y(kkk),'*')
for i = 1:numel(num)
    name = strcat(prefix_name,middle_name,num2str(num(i),'%04d'),'.h5');
    [data,domain] = getData(char(name),char('/Flbm_0'),2,true);
    
    Flbmx = reshape(data(:,:,:,1),[Nx,Ny]);
    Flbmy = reshape(data(:,:,:,2),[Nx,Ny]);
    Flbm(i,1) = sum(sum(Flbmy));
    Flbmxx(i,1) = sum(sum(Flbmx));
    [data,domain] = getData(char(name),char('/Velocity_0'),2,true);
    vx = reshape(data(:,:,:,1),[Nx,Ny]);
    vy = reshape(data(:,:,:,2),[Nx,Ny]);
    v = (vx.^2 + vy.^2).^0.5;
    U(i,1) = mean(v(kkk));
    [data,domain] = getData(char(name),char('/Density_0'),2,false);
    rho = reshape(data,[Nx,Ny]);
    Rho(i,1) = mean(rho(kkk));
end
Cd = (Flbm./((U.^2).*Rho))./((0.5).*2*R)
plot(-Cd,'-*')