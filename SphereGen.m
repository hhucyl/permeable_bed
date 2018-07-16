clear
clc
prefix_name = {'/home/pzhang/chen/premeable_bed/'};
file_name = {'sphere.h5'};
start = 16;
R = 10;
xs = 50;ys=50;zs=50;
l = 100;
%for xs = 20:10:70
    ys = xs;
    zs = xs;
    xx = [xs xs+l];
    yy = [ys ys+l];
    zz = [zs zs+l];
    pos = getData(char(file_name),char('/Position'),1,false);
    ppos = pos(start:end);
    num = numel(ppos)/3;
    px = ppos(1:3:end-2);
    py = ppos(2:3:end-1);
    pz = ppos(3:3:end);
    kkkx = find(px>=xx(1) & px<=xx(2));
    kkky = find(py>=yy(1) & py<=yy(2));
    kkkz = find(pz>=zz(1) & pz<=zz(2));
    kkk = intersect(kkkx,kkky);
    kkk = intersect(kkk,kkkz);
    temp = [px(kkk)-xx(1),py(kkk)-yy(1),pz(kkk)-zz(1)];
    figure
    sphereplot(temp(:,1),temp(:,2),temp(:,3),R);
    n = numel(kkk)*4.0*3.0*pi*R^3/l^3;
%     fid = fopen('sphere.txt','w');
%     fprintf(fid,'%d\n',numel(kkk));
%     fprintf(fid,'%f %f %f\n',temp);
%     fclose(fid)
    A = 0;
%     for i = 1:numel(kkk)
%         flagx = temp(i,1)-R>0 & temp(i,1)+R<l;
%         flagy = temp(i,2)-R>0 & temp(i,2)+R<l;
%         flagz = temp(i,3)-R>0 & temp(i,3)+R<l;
%         if(flagx&flagy&flagz)
%             A = A + 4.0*pi*R^2;
%         end
%     end
%end