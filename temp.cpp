#include "./lbm/Domain.h"


void Initial(LBM::Domain &dom, double rho, Vec3_t &v0, Vec3_t&g0)
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
        
    }
    dom.Rho0 = rho;
}



using namespace std;

int main (int argc, char **argv) try
{
      
    size_t h = 20;
    size_t nx = h;
    size_t ny = h;
    size_t nz = 1;//0.12/0.038 = 3.1579
    double dx = 1.0;
    double dt = 1.0;
    double nu = 0.1;
    Vec3_t pos(nx/2-1,ny/2-1,0.0);
    double XX = 0.85;//x = (c/cmax)^1/3
    double cmax = M_PI/6.0;//for bcc cmax = 3^0.5*pi/8
    double R = std::pow(XX*XX*XX*cmax*h*h*h/(4.0/3.0*M_PI),1.0/3.0);//x = (c/cmax)^1/3 c = Vs/V
    std::cout<<"R = "<<R<<std::endl;
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(0.0,0.0,2e-5);
    std::cout<<"nu "<<nu<<std::endl;
    std::cout<<"Nx "<<nx<<std::endl;
    std::cout<<"Ny "<<ny<<std::endl;
    std::cout<<"Nz "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT,nu,iVec3_t(nx,ny,1),dx,dt);
    dom.Nproc = 1;
    dom.Sc = 0.0;


    //initial
    Initial(dom,rho0,v0,g0);
   
    double Tf = 2;
    double dtout = 1;
    char const * TheFileKey = "temp";
    //solving
    dom.StartTime = std::clock();
    dom.StartSolve();
    double tout = 0;
    while(dom.Time<Tf)
    {
        
        if (dom.Time>=tout)
        {
            
            String fn;;
            fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            dom.WriteXDMF(fn.CStr());
            dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" ---"<<std::endl;
            tout += dtout;
        }
        // (dom.*dom.ptr2collide)();
        
        dom.SetZero();
        dom.AddDiskG(pos,5.0);
        dom.CollideMRT();
        // InOut(dom,&my_dat);
        dom.Stream();
        
        dom.CalcProps();
        dom.Time += 1;
    }
    // dom.EndSolve();
    
    dom.EndTime = std::clock();
    double ttime = (double)(dom.EndTime - dom.StartTime)/CLOCKS_PER_SEC;
    printf("\033[01;34m---Elapsed Time = %f s---\033[0m\n",ttime);
    
    return 0;
}MECHSYS_CATCH