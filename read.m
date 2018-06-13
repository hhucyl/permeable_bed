clear
clc
pname1 = {'test_carve3b_'};
num = 1:2;
for i = 1:numel(num)
   
    
    name = strcat(pname1,num2str(num(i),'%04d'),'.h5');
    [data, domain] = getData(char(name),char('/F_0'),2,false);
    Fb(:,i) = data;
    Dim(1) = domain.Nx;
    Dim(2) = domain.Ny;
    Dim(3) = domain.Nz;
    [Fbt(:,i),domain] = getData(char(name),char('/Ftemp_0'),2,false);
    [data, domain] = getData(char(name),char('/Gamma'),2,false);
    gama = data;
    [Velb(:,i),domain] = getData(char(name),char('/Velocity_0'),2,false);
end
kkk = find(gama<2);
diff = [];
diff = Fb(:,1) - Fbt(:,1);
C = {
	{ 0, 0, 0},... 
    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},... 
    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 0, 1}, {-1, 0,-1},...
    { 1, 0,-1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1}
};
for i= 1:(Dim(1)*Dim(2)*Dim(3))
% i=1;
    x1 = (i-1)*19 + 1;
    x2 = 19*i;
    f = Fbt(x1:x2,1);
    ttemp = [0 0 0];
    for j = 1:19
        c{j,1} = [C{j}{1},C{j}{2},C{j}{3}];
        ttemp = ttemp + f(j).*c{j,1};
    end
    v{i,1} = ttemp;
    x1 = (i-1)*3 + 1;
    x2 = i*3;
    vv{i,1} = Velb(x1:x2,1);
end

i=50
ix = mod(i,Dim(1));
iy = mod(floor(i/Dim(1)),Dim(2));
iz = floor(i/(Dim(1)*Dim(2)));