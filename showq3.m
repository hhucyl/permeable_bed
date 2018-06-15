clear
clc
prefix_name = {'/home/pzhang/chen/permeable_bed/'};
middle_name = {'test_carve3b_'};

name = strcat(prefix_name,middle_name,num2str(0,'%04d'),'.h5');
[data,domain] = getData(char(name),char('/q_0'),2,false);
Nneigh = 19;
Nx = domain.Nx;
Ny = domain.Ny;
Nz = domain.Nz;
count = 1;
% for k = 1:Nx
% for j = 1:Ny
for i = 1:Nneigh
    q{i} = data(i:Nneigh:(end-Nneigh+i));
    q{i} = reshape(q{i},[Nx,Ny,Nz]);
%       is = 1+(count-1)*Nneigh;
%       ie = Nneigh*count;
%       q{i,j,k} = data(is:ie);
%       count = count+1;
end
% end
% end
[datag,domain] = getData(char(name),char('/Gamma'),2,false);
gamma = reshape(datag,[Nx,Ny,Nz]);
if(Nneigh == 15)
    C = {[ 0, 0, 0];
         [1, 0, 0];
         [-1, 0, 0]; 
         [ 0, 1, 0]; 
         [ 0,-1, 0];
         [ 0, 0, 1]; 
         [ 0, 0,-1]; 
         [ 1, 1, 1];
         [-1,-1,-1]; 
         [ 1, 1,-1]; 
         [-1,-1, 1];
         [ 1,-1, 1]; 
         [-1, 1,-1]; 
         [ 1,-1,-1];
         [-1, 1, 1]};
end
if(Nneigh == 19)
    C = {[ 0, 0, 0]; 
         [ 1, 0, 0]; 
         [-1, 0, 0]; 
         [ 0, 1, 0]; 
         [ 0,-1, 0]; 
         [ 0, 0, 1]; 
         [ 0, 0,-1]; 
         [ 1, 1, 0]; 
         [-1,-1, 0]; 
         [ 1, 0, 1]; 
         [-1, 0,-1]; 
         [ 0, 1, 1]; 
         [ 0,-1,-1];
         [ 1,-1, 0]; 
         [-1, 1, 0]; 
         [ 1, 0,-1]; 
         [-1, 0, 1]; 
         [ 0, 1,-1]; 
         [ 0,-1, 1]};
end
    
pos = [25.5,105.5,25.5];
R = 5;
ipt = [];
for k = 2:Nneigh
%     k=6;
    temp = q{k};
    I = find(temp>0);
    II{k-1} = I;
    [ix iy iz] = ind2sub(size(temp),I);
    c = C{k};
    ix = ix-1;
    iy = iy-1;
    iz = iz-1;
    ixx = ix-pos(1);
    iyy = iy - pos(2);
    izz = iz - pos(3);
    nix = ix+c(1);
    niy = iy+c(2);
    niz = iz+c(3);
    nixx = nix-pos(1);
    niyy = niy-pos(2);
    nizz = niz-pos(3);
    ipx = ix+temp(I).*c(1);
    ipy = iy+temp(I).*c(2);
    ipz = iz+temp(I).*c(3);
    ipxx = ipx -pos(1);
    ipyy = ipy - pos(2);
    ipzz = ipz - pos(3);
    diff(k,1) = sum(ipxx.^2+ipyy.^2+ipzz.^2 - R*R);
%     scatter3(ixx,iyy,izz,'ro')
%     hold on
%     scatter3(nixx,niyy,nizz,'ro')
    scatter3(ipxx,ipyy,ipzz,'b*')
    ipt = [ipt;ipxx,ipyy,ipzz];
    axis equal
    drawnow
    hold on
    for i = 1:numel(I)
        xx = [ixx(i),nixx(i)];
        yy = [iyy(i),niyy(i)];
        zz = [izz(i),nizz(i)];
        plot3(xx,yy,zz,'.r-')
        hold on
%         scatter3(ix(i),iy(i),iz(i),'bs')
%         hold on
% %         axis equal
% %         drawnow
    end
%     if(qq>0)
%         if s<2 & s>0
%             plot3(ix,iy,iz,'bs')
%         end
%         hold on
%         c = C{k};
%         ipx = ix+qq*c(1);
%         ipy = iy+qq*c(2);
%         ipy = iz+qq*c(3);
%         plot3(ipx,ipy,ipz,'*')
%         s = s+1;
%     end
end
plot3(0,0,0,'kd')
         
%axis([180 220 250 350])
 
 