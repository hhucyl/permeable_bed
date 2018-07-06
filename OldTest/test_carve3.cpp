#include "./lbm/Domain.h"


struct myUserData
{
    double u0;
    double Tff;
};

void Setup(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    size_t Nproc = dom.Nproc;
    double k = 1.0;
    // if(dom.Time<=dat.Tff)
    // {
    //     k = dom.Time/dat.Tff;
    // }
    size_t index = 0;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t iy=0; iy<ny; ++iy)
    for(size_t iz=0; iz<nz; ++iz)
    {
    // for(size_t i=0; i<dom.Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     dom.idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        if(dom.IsSolid[index][iy][iz]) continue;            
        // double ux = dat.u0*k;
        // // double uy = -dat.u0*(iz)*(nz-1-iz)/(0.5*nz*nz);
        // double uy = 0.0;
        // double uz = 0.0;
        // double *f = dom.F[index][iy][iz];
        
        

        // double *f1 = dom.F[index+1][iy][iz];
        // double rho1 = dom.Rho[index+1][iy][iz];
        // Vec3_t vel(ux,uy,uz);
        // Vec3_t vel1 = dom.Vel[index+1][iy][iz];
        // for(size_t k=0; k<dom.Nneigh; ++k)
        // {
        //     f[k] = dom.Feq(k,rho1,vel)+f1[k] - dom.Feq(k,rho1,vel1);
        // }
        double *f = dom.F[index][iy][iz];
        Vec3_t vv(dat.u0,0.0,0.0);
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,1.0, vv);
        }

    }
    
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t iy=0; iy<ny; ++iy)
    for(size_t iz=0; iz<nz; ++iz)
    { 
    // for(size_t i=0; i<dom.Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     dom.idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        index = nx-1;
        if(dom.IsSolid[index][iy][iz]) continue; 
        double *f = dom.F[index][iy][iz];
        double *f1 = dom.F[index-1][iy][iz];
        
        // dom.Vel[ix][index][iz] = dom.Vel[ix][index+1][iz];
        //dom.Rho[ix][index][iz] = dom.Rho[ix][index+1][iz];
        for (size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = f1[k];
        }
        
    // }
    // for(size_t i=0; i<dom.Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     dom.idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        // index = 0;
        // if(dom.IsSolid[ix][index][iz]) continue; 
        // double *f = dom.F[ix][index][iz];
        // double *f1 = dom.F[ix][index+1][iz];

        // double rho = dom.Rho[ix][index+1][iz];
        // Vec3_t vel1 = dom.Vel[ix][index+1][iz];
        // for (size_t k=0; k<dom.Nneigh; ++k)
        // {
        //     f[k] = dom.Feq(k,1.0,vel1)+f1[k] - dom.Feq(k,rho,vel1);
        // }
        
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
    for(size_t ix=0; ix<dom.Ndim(0); ix++)
    for(size_t iy=0; iy<dom.Ndim(1); iy++)
    for(size_t iz=0; iz<dom.Ndim(2); iz++)
    {
        dom.Rho[ix][iy][iz] = rho;
        dom.Vel[ix][iy][iz] = 0.0, 0.0, 0.0;
        dom.BForce[ix][iy][iz] = 0.0,1e-5,0.0;
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
    
    size_t Nproc = 12;
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 100;
    double dx = 1.0;
    double dt = 1.0;
    double u0 = 0.05;
    double Re = 30;
    double R = 5;
    double nu = u0*R/Re;
    nu = 0.1;
    Vec3_t pos(49.5,49.5,49.5);
    // Vec3_t pos(200.5,299.5,0);
    
    
    LBM::Domain dom(D3Q19,MRT,LIBB,nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.u0 = u0;
    dom.Nproc = Nproc;

    dom.AddSphereQ(pos,R);
    dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    //bounndary
    // dom.AddDisk(pos,R);
    //initial
    double rho = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0);
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
    double t = 500;
    dom.Solve(Tf,t,NULL,NULL,"test_carve3b",true,Nproc);
    return 0;
}MECHSYS_CATCH