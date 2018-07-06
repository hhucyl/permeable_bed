#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
};

void Report(LBM::Domain &fluid_dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = fluid_dom.Ndim(0);
    size_t ny = fluid_dom.Ndim(1);
    size_t nz = fluid_dom.Ndim(2);
    double dx = fluid_dom.dx;
    if(fluid_dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","permeability");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"P1"<<Util::_8s<<"P2"<<Util::_8s<<"U"<<Util::_8s<<"K\n";
    }else{
        size_t index = 0;
        double P1 = 0;
        double num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        {
            // if(fluid_dom.IsSolid[ix][iy][index]) continue;
            // num += 1;
            P1 += fluid_dom.Rho[ix][iy][index]/3.0;
        }
        P1 /= nx*ny;
        index = nz -1;

        double P2 = 0;
        num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        {
            // if(fluid_dom.IsSolid[ix][iy][index]) continue;
            // num +=1.0;
            P2 += fluid_dom.Rho[ix][iy][index]/3.0;
        }
        P2 /= nx*ny;
        double U = 0;
        num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        for(size_t iz=0;iz<nz;++iz)
        {
            // if(fluid_dom.IsSolid[ix][iy][iz]) continue;
            // num += 1.0;
            Vec3_t e(0,0,1);
            U += dot(fluid_dom.Vel[ix][iy][0],e);
        }
        U /= nx*ny*nz;
        // U /= num;
        U += 0.5*dat.g;
        // std::cout<<dat.g<<" "<<dat.R<<std::endl;
        // double K = -(U*dat.nu)/((P2-P1)/((double)nz) - dat.g);
        double K = (U*dat.nu)/(dat.g);
        // double K = 8.0*dat.nu*nz*nz*U/(9*nz*(P1-P2));
        K = (6.0*M_PI*dat.R*K)/(nx*ny*nz);
        dat.oss_ss<<Util::_10_6<<fluid_dom.Time<<Util::_8s<<P1<<Util::_8s<<P2<<Util::_8s<<U<<Util::_8s<<K<<std::endl;
    }
}
void Initial(LBM::Domain &dom, double rho, Vec3_t &v0,  Vec3_t &g0)
{
    for(size_t ix=0; ix<dom.Ndim(0); ix++)
    for(size_t iy=0; iy<dom.Ndim(1); iy++)
    for(size_t iz=0; iz<dom.Ndim(2); iz++)
    {
        dom.Rho[ix][iy][iz] = rho;
        dom.Vel[ix][iy][iz] = 0.0, 0.0, 0.0;
        dom.BForce[ix][iy][iz] = g0;
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.F[ix][iy][iz][k] = dom.Feq(k,rho,v0);            
            dom.Ftemp[ix][iy][iz][k] = dom.Feq(k,rho,v0);            
        }
    // std::cout<<dom.F[ix][iy][iz][18]<<std::endl;
        
    }
    // dom.Rho0 = rho;//very important
}

int main () try
{
    
    size_t Nproc = 8;
    size_t h = 30;
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double dx = 1.0;
    double dt = 1.0;
    double XX = 0.85;//x = (c/cmax)^1/3
    double cmax = M_PI/6.0;//for bcc cmax = 3^0.5*pi/8
    double R = std::pow(XX*XX*XX*cmax*nx*nx*nx/(4.0/3.0*M_PI),1.0/3.0);//x = (c/cmax)^1/3 c = Vs/V
    std::cout<<"R = "<<R<<std::endl;
    double ddx = -0.5;
    Vec3_t pos(ddx,ddx,ddx);
    // Vec3_t pos(200.5,299.5,0);
    double nu = (0.6-0.5)/3.0;
    //nu = 1.0/30.0;
    LBM::Domain dom(D3Q19,MRT,LIBB,nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = 2e-5;
    my_dat.R = R;
    Vec3_t g0(0.0,0.0,my_dat.g);
    dom.Nproc = Nproc;


    dom.AddSphereQ(pos,R);
    pos = nx-1-ddx,0+ddx,0+ddx;
    dom.AddSphereQ(pos,R);
    pos = 0+ddx,ny-1-ddx,0+ddx;
    dom.AddSphereQ(pos,R);
    pos = nx-1-ddx, ny-1-ddx, 0+ddx;
    dom.AddSphereQ(pos,R);
    //pos = nx/2.0-1,ny/2.0-1,nz/2.0-1;
    //dom.AddSphereQ(pos,R);
    pos = 0.0+ddx,0.0+ddx,nz-1-ddx;
    dom.AddSphereQ(pos,R);
    pos = nx-1-ddx,0+ddx,nz-1-ddx;
    dom.AddSphereQ(pos,R);
    pos = 0+ddx,ny-1-ddx,nz-1-ddx;
    dom.AddSphereQ(pos,R);
    pos = nx-1-ddx, ny-1-ddx, nz-1-ddx;
    dom.AddSphereQ(pos,R);             

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    //bounndary
    // dom.AddDisk(pos,R);
    //initial
    double rho = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0,g0);
    // for(size_t ix=0; ix<nx; ++ix)
    // for(size_t iy=0; iy<ny; ++iy)
    // {
    //     dom.IsSolid[ix][iy][0] = true;
    //     dom.IsSolid[ix][iy][nz-1] = true;
    // }
    /*for(size_t iz=0; iz<nz; ++iz)
    for(size_t iy=0; iy<ny; ++iy)
    {
        dom.IsSolid[0][iy][iz] = true;
        dom.IsSolid[nx-1][iy][iz] = true;
    }*/

    

    double Tf = 100000;
    // my_dat.Tff = 0.7*Tf;
    double t = 100;
    dom.Solve(Tf,t,NULL,Report,"test_bcc",true,Nproc);
    return 0;
}MECHSYS_CATCH