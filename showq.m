clear
clc
prefix_name = {'/home/pzhang/chen/permeable_bed/'};
middle_name = {'test_carve_'};

name = strcat(prefix_name,middle_name,num2str(1,'%04d'),'.h5');
[data,domain] = getData(char(name),char('/q_0'),2,false);
Nneigh = 9;
Nx = domain.Nx;
Ny = domain.Ny;
Nz = domain.Nz;
for i = 1:Nneigh
    q{i} = data(i:Nneigh:(end-Nneigh+i));
    q{i} = reshape(q{i},[Nx,Ny,Nz]);
end
C = {[0,0,0];
     [1,0,0];
     [0,1,0];
     [-1,0,0];
     [0,-1,0];
     [1,1,0];
     [-1,1,0];
     [-1,-1,0];
     [1,-1,0]};
 for ix = 180:220
     for iy = 250:350
         
         s = 0;
         for k = 1:Nneigh
             temp = q{k};
             qq = temp(ix,iy);
             if(qq>0)
                if s<2 & s>0
                    plot(ix,iy,'bs')
                end
                hold on
                c = C{k};
                ipx = ix+qq*c(1);
                ipy = iy+qq*c(2);
                plot(ipx,ipy,'*')
                s = s+1;
             end
             
         end
         axis([180 220 250 350])
         %axis equal
         drawnow
     end
 end