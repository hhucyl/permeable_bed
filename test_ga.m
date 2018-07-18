clear
clc
prefix_name = {'/home/pzhang/chen/permeable_bed/'};
middle_name = {'test_line_'};
num = 1:300;
q = 0.25;
if q<0.5
    ga = 0.5 -q;
else
    ga = 1.5 -q;
end
tau = 0.8;
nu = (tau-0.5)/3.0;
for i = 1:numel(num)
%     i = 2
    name = strcat(prefix_name,middle_name,num2str(num(i),'%04d'),'.h5');
    [data, domain] =getData(char(name),char('/Velocity_0'),2,true);
    Nx = domain.Nx;
    Ny = domain.Ny;
    Nz = domain.Nz;
    u = data(:,:,Nz/2,1);
    data = getData(char(name),char('/BForce_0'),2,false);
    g = data(1);
    data = getData(char(name),char('/Velocity_P_0'),2,false);
    uw = data(end-2);
    
    H = Ny-1+q-0.5;
    yr = 0:0.01:H;
    y = 0:1:Ny-1;
    y = y-0.5;
    urr = -g/(2.0*nu).*y.^2 + (g/(2.0*nu)*H+uw/H).*y;
    ur = -g/(2.0*nu).*yr.^2 + (g/(2.0*nu)*H+uw/H).*yr;
    plot(ur,yr,'r')
    hold on
    uu =  u(Nx/2,:);
    plot(uu,y,'o')
    hold off
    drawnow
    r = (uu-urr)./urr;
    rr(i,1) = sum(r(2:end));
end
figure
plot(rr,'-*')