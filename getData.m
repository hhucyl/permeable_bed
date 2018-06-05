function [data, domain] = getData(file_name,data_name,type,flag3d)
    if(type == 1)
        %domain.Np = double(h5read(char(file_name)),'/NP');
        if(flag3d)
            datat = h5read(char(file_name),char(data_name));
            data(:,1) = datat(1:3:(end-2));
            data(:,2) = datat(2:3:(end-1));
            data(:,3) = datat(3:3:(end));
        else
            data = h5read(char(file_name),char(data_name));
        end
    end
    if(type == 2)
        domain.Nx = double(h5read(char(file_name),'/Nx'));
        domain.Ny = double(h5read(char(file_name),'/Ny'));
        domain.Nz = double(h5read(char(file_name),'/Nz'));
        Nx = domain.Nx;
        Ny = domain.Ny;
        Nz = domain.Nz;
        if(flag3d)
                datat = h5read(char(file_name),char(data_name));
                Vx = datat(1:3:(end-2));
                Vy = datat(2:3:(end-1));
                Vz = datat(3:3:end);
                vx = reshape(Vx,[Nx,Ny,Nz]);
                vy = reshape(Vy,[Nx,Ny,Nz]);
                vz = reshape(Vz,[Nx,Ny,Nz]);
                data(:,:,:,1) = vx;
                data(:,:,:,2) = vy;
                data(:,:,:,3) = vz;
        else
            data = h5read(char(file_name),char(data_name));
%             data = reshape(data,[domain.Nx,domain.Ny,domain.Nz]);
        end
    end
        
end