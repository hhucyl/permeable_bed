#include "./lbm/Domain.h"


struct myUserData
{
    double u0;
};

void Setup(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.nx;
    size_t ny = dom.ny;
    size_t nz = dom.nz;
    size_t Nproc = dom.Nproc;
    size_t index = ny-1;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iz=0; iz<nz; ++iz)
    {   
        double ux = 0.0;
        double uy = -dat.u0;
        double *f = dom.F[ix][index][iz];
        
        double rho = 1.0/(1.0-uy)*(f[0]+f[1]+f[3]+2.0*(f[2]+f[5]+f[6]));
        dom.Vel[ix][index][iz] = ux,uy,0.0;
        dom.Rho[ix][index][iz] = rho;
        f[4] = f[2]- 2.0/3.0 * dom.Rho[ix][index][iz]*uy;
        f[8] = f[6] + 0.5*(f[3] - f[1]) + 0.5*rho*ux - 1.0/6.0*rho*uy;
        f[7] = f[5] + 0.5*(f[1] - f[3]) - 0.5*rho*ux - 1.0/6.0*rho*uy;

        // double *f1 = dom.F[ix][index-1][iz];
        // double rho = dom.Rho[ix][index-1][iz];
        // Vec3_t vel(ux,uy,0.0);
        // Vec3_t vel1 = dom.Vel[ix][index-1][iz];
        // for(size_t k=0; k<dom.Nneigh; ++k)
        // {
        //     f[k] = dom.Feq(k,rho,vel)+f1[k] - dom.Feq(k,rho,vel1);
        // }

    }
    
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iz=0; iz<nz; ++iz)
    {  
        index = 0;
        double ux = 0.0;
        double uy = -dat.u0;
        double *f = dom.F[ix][index][iz];
        double *f1 = dom.F[ix][index+1][iz];
        
        // dom.Vel[ix][index][iz] = dom.Vel[ix][index+1][iz];
        //dom.Rho[ix][index][iz] = dom.Rho[ix][index+1][iz];
        for (size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = f1[k];
        }
        
    }
    // #ifdef USE_OMP
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    // #endif
    // for(size_t ix=0; ix<nx; ++ix)
    // for(size_t iz=0; iz<nz; ++iz)
    // {  
    //     index = 0;
    //     double ux = 0.0;
    //     double uy = -dat.u0;
    //     double *f = dom.F[ix][index][iz];
        
    //     double rho = 1.0/(1.0-uy)*(f[0]+f[1]+f[3]+2.0*(f[4]+f[8]+f[7]));
    //     dom.Vel[ix][index][iz] = ux,uy,0.0;
    //     dom.Rho[ix][index][iz] = rho;
    //     f[2] = f[4] + 2.0/3.0 * dom.Rho[ix][index][iz]*uy;
    //     f[6] = f[8] - 0.5*(f[3] - f[1]) - 0.5*rho*ux + 1.0/6.0*rho*uy;
    //     f[5] = f[7] - 0.5*(f[1] - f[3]) + 0.5*rho*ux + 1.0/6.0*rho*uy;
    // }
}

void Initial(LBM::Domain &dom, double rho, Vec3_t &v0)
{
    for(size_t ix=0; ix<dom.nx; ix++)
    for(size_t iy=0; iy<dom.ny; iy++)
    for(size_t iz=0; iz<dom.nz; iz++)
    {
        dom.Rho[ix][iy][iz] = rho;
        dom.Vel[ix][iy][iz] = v0;
        dom.Check[ix][iy][iz] = 0.0;
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.F[ix][iy][iz][k] = dom.Feq(k,rho,v0);            
        }
    }
}

int main () try
{
    LBM::Domain dom;
    size_t Nproc = 10;
    size_t nx = 400;
    size_t ny = 400;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double u0 = 0.05;
    double Re = 100;
    double R = 20;
    double nu = u0*R/Re;
    Vec3_t pos(200.5,299.5,0);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.u0 = u0;
    dom.AllocateDomainMem(nx,ny,1,nu,dx,dt);
    dom.Nproc = Nproc;
    dom.Isq = true;
    //bounndary
    dom.AddDisk(pos,R);
    //initial
    double rho = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0);
    
    
    double Tf = 1e6;
    double t = Tf/1e3;
    dom.Solve(Tf,t,Setup,NULL,"test_carve",true);
    return 0;
}MECHSYS_CATCH