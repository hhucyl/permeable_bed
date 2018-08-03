#include "./lbm/Domain.h"

struct myUserData
{
    std::ofstream oss_ss;
    double R;
    double u;
};

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

void InOut(LBM::Domain &dom,  void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    size_t Nneigh = dom.Nneigh;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
    for(size_t iz=0; iz<nz; iz++)
    {
        for(size_t iy=0; iy<ny; iy++)
        {
            double *f = dom.Ftemp[0][iy][iz];
            double *f1 = dom.Ftemp[1][iy][iz];
            double rho1 = dom.Rho[1][iy][iz];
            Vec3_t vel1 = dom.Vel[1][iy][iz];
            double rho = dom.Rho[0][iy][iz];
            Vec3_t vel = dom.Vel[0][iy][iz];
            
            if(dom.IsSolid[0][iy][iz]) continue;
            
            double uu = dat.u;
            Vec3_t vv (uu,0.0,0.0);

            for(size_t k=0; k<Nneigh; k++)
            {
                // if(k==1&&k==5&&k==7&&k==13&&k==15)
                f[k] = dom.Feq(k,rho1,vv) + f1[k] - dom.Feq(k,rho1,vel1);
                // f[k] = dom.Feq(k,rho1,vv);
            }
            //BForce[ix][iy][iz] = OrthoSys::O;
            vel = OrthoSys::O;
            rho = 0.0;
            for (size_t k=0;k<Nneigh;k++)
            {
                rho +=  f[k];
                vel +=  f[k]*dom.C[k];
            }
            vel *= dom.Cs/rho;

        }
    }
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
    for(size_t iz=0; iz<nz; iz++)
    {
        for(size_t iy=0; iy<ny; iy++)
        {
            double *f = dom.Ftemp[nx-1][iy][iz];
            double *f1 = dom.Ftemp[nx-2][iy][iz];
            double rho = dom.Rho[nx-1][iy][iz];
            Vec3_t vel = dom.Vel[nx-1][iy][iz];
            if(dom.IsSolid[nx-1][iy][iz]) continue;
           

            for(size_t k=0; k<Nneigh; k++)
            {
                f[k] = f1[k];
            }
            //BForce[ix][iy][iz] = OrthoSys::O;
            vel = OrthoSys::O;
            rho = 0.0;
            for (size_t k=0;k<Nneigh;k++)
            {
                rho +=  f[k];
                vel +=  f[k]*dom.C[k];
            }
            vel *= dom.Cs/rho;
            

        }
    }

}

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    double nu = dom.Nu;
    if(dom.Time <1e-6)
    {
        String fs;
        // fs.Printf("%s_%d_%g.out","R",dat.bbtype,dom.Tau);
        fs.Printf("%s.out","Cd");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"U"<<Util::_8s<<"Re"<<Util::_8s<<"Cd\n";
    }else{
        double F = 0.0;
        double rho = 0.0;
        for(size_t iz=0; iz<nz; ++iz)
        for(size_t ix=0; ix<nx; ++ix)
        for(size_t iy=0; iy<ny; ++iy)
        {
            F += dom.Flbm[ix][iy][iz](0);
            rho += dom.Rho[ix][iy][iz];
        }
        rho /= nx*ny*nz;
        double Re = dat.u*2*dat.R/nu;
        double Cd = 2.0*F/(rho*dat.u*dat.u*dat.R*dat.R*M_PI);
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dat.u<<Util::_8s<<Re<<Util::_8s<<Cd<<std::endl;
        
    }
}

using namespace std;

int main (int argc, char **argv) try
{
      
    size_t Nproc = 12;
    double Tf = 1e5; 
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 100;//0.12/0.038 = 3.1579
    double dx = 1.0;
    double dt = 1.0;
    double nu = 0.01;
    Vec3_t pos(nx/2-1,ny/2-1,nz/2-1);
    double R = 5;
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(0.0,0.0,0.0);
    std::cout<<"nu "<<nu<<std::endl;
    std::cout<<"Nx "<<nx<<std::endl;
    std::cout<<"Ny "<<ny<<std::endl;
    std::cout<<"Nz "<<nz<<std::endl;
    LBM::Domain dom(D3Q19,SRT,nu,iVec3_t(nx,ny,nz),dx,dt);
    dom.Nproc = Nproc;
    dom.Sc = 0.0;
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.R = R;
    my_dat.u = 0.04;
    //initial
    Initial(dom,rho0,v0,g0);
    std::fstream ifile("Nodes.txt",std::ios::in);
	int N = 0;
    double Rc = 0;
    Vec3_t post(0.0,0.0,0.0);
	if(!ifile.fail())
	{
		ifile>>Rc;
        if(std::fabs(Rc-R)>1e-6)
        {
            throw new Fatal("Radius is NOT RIGHT"); 
        }
        ifile>>N;
        double ds = (M_PI*dx)/(3.0*N)*(12.0*Rc*Rc + dx*dx);
		for(size_t i=0;i<N;++i)
		{
			double x;
			double y;
			double z;
    		ifile>>x>>y>>z;
            post = x+pos(0), y+pos(1), z+pos(2);
            if(std::fabs(dot(post - pos, post - pos)-25)>1e-5) 
            {
                std::cout<<std::fabs(dot(post - pos, post - pos)-25)<<std::endl;
            }
            dom.points.push_back(post);
            dom.dS.push_back(ds);
		}
	} 

    Tf = 1e5;
    double dtout = 1e2;
    char const * TheFileKey = "test_ibm";
    //solving
    dom.StartTime = std::clock();
    dom.StartSolve();
    double tout = 0;
    while(dom.Time<Tf)
    {
        
        if (dom.Time>=tout)
        {
            
            String fn;
            // fn.Printf("%s_%d_%g_%04d", TheFileKey,bbtype, tau,dom.idx_out);
            fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            dom.WriteXDMF(fn.CStr());
            dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        // (dom.*dom.ptr2collide)();
        
        
        dom.ApplyIBM3D(pos,R);
        dom.CollideSRTIBM();
        InOut(dom,&my_dat);
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