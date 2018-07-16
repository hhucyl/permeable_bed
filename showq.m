clear
clc
prefix_name = {'/home/pzhang/chen/permeable_bed/'};
middle_name = {'test_line_'};

name = strcat(prefix_name,middle_name,num2str(1,'%04d'),'.h5');
[data,domain] = getData(char(name),char('/q_0'),2,false);
Nneigh = 19;
Nx = domain.Nx;
Ny = domain.Ny;
Nz = domain.Nz;
for i = 1:Nneigh
    q{i} = data(i:Nneigh:(end-Nneigh+i));
    q{i} = reshape(q{i},[Nx,Ny,Nz]);
end
if(Nneigh == 9)
    C = {[0,0,0];
         [1,0,0];
         [0,1,0];
         [-1,0,0];
         [0,-1,0];
         [1,1,0];
         [-1,1,0];
         [-1,-1,0];
         [1,-1,0]};
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
 ipt = [];
for k = 2:Nneigh
%      k=4;
    temp = q{k};
    I = find(temp>0);
    II{k-1} = I;
    [ix iy iz] = ind2sub(size(temp),I);
    c = C{k};
    ix = ix-1;
    iy = iy-1;
    iz = iz-1;
    ixx = ix;
    iyy = iy ;
    izz = iz ;
    nix = ix+c(1);
    niy = iy+c(2);
    niz = iz+c(3);
    nixx = nix;
    niyy = niy;
    nizz = niz;
    ipx = ix+temp(I).*c(1);
    ipy = iy+temp(I).*c(2);
    ipz = iz+temp(I).*c(3);
%     scatter3(ixx,iyy,izz,'ro')
%     hold on
%     scatter3(nixx,niyy,nizz,'ro')
    if(Nneigh == 9)
        plot(ipx,ipy,'b*')
    end
    if(Nneigh == 19)
        plot3(ipx,ipy,ipz,'b*')
    end
    ipt = [ipt;ipx,ipy,ipz];
    axis equal
    drawnow
    hold on
    for i = 1:numel(I)
        xx = [ixx(i),nixx(i)];
        yy = [iyy(i),niyy(i)];
        zz = [izz(i),nizz(i)];
        if(Nneigh == 9)
            plot(xx,yy,'.r-')
        end
        if(Nneigh == 19)
            plot3(xx,yy,zz,'.r-')
        end
        if(Nneigh == 9)
            plot(xx(2),yy(2),'.k','markersize',20)
        end
        if(Nneigh == 19)
            plot3(xx(2),yy(2),zz(2),'.k','markersize',20)
        end
        hold on
%         scatter3(ix(i),iy(i),iz(i),'bs')
%         hold on
% %         axis equal
        drawnow
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