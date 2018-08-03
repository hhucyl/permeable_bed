#include "./lbm/Domain.h"

struct myUserData
{
    std::ofstream oss_ss;
    double R;
    double u;
    double g;
    double nu;
    double Tf;
    double K;
    int bbtype;
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
    double dx = dom.dx;
    if(dom.Time <1e-6)
    {
        String fs;
        // fs.Printf("%s.out","permeability");
        fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"U"<<Util::_8s<<"r"<<Util::_8s<<"K\n";
    }else{
        size_t index = 0;
        double U = 0;
        double num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        for(size_t iz=0;iz<nz;++iz)
        {
            // if(dom.IsSolid[ix][iy][iz]) continue;
            // num += 1.0;
            Vec3_t e(0,0,1);
            U += dot(dom.Vel[ix][iy][iz],e);
        }
        U /= nx*ny*nz;
        // U /= num;
        // U += 0.5*dat.g;
        // std::cout<<dat.g<<" "<<dat.R<<std::endl;
        // double K = -(U*dat.nu)/((P2-P1)/((double)nz) - dat.g);
        double K = (U*dat.nu)/(dat.g);
        // double K = 8.0*dat.nu*nz*nz*U/(9*nz*(P1-P2));
        K = (6.0*M_PI*dat.R*K)/(nx*ny*nz);
        double r = std::fabs(K-dat.K)/dat.K;
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<U<<Util::_8s<<r<<Util::_8s<<K<<std::endl;
        
        if(r<1e-5)
        {
            dom.Time = dat.Tf;
        }
        dat.K = K;
    }
}

using namespace std;

int main (int argc, char **argv) try
{
      
    size_t Nproc = 12;
    double Tf = 1e5;
    size_t h = 30;
    size_t h1 = 30;
    int bbtype = -2;
    double tau = 0.6;
    if(argc>=2) bbtype = atoi(argv[1]); 
    if(argc>=3) h = atoi(argv[2]);
    if(argc>=4) tau = atof(argv[3]);     
    if(argc>=5) Nproc = atoi(argv[4]);  
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;//0.12/0.038 = 3.1579
    double dx = 1.0;
    double dt = 1.0;
    double nu = 0.1;
    nu = (tau-0.5)/3.0;
    Vec3_t pos(nx/2-1,ny/2-1,nz/2-1);
    double XX = 0.85;//x = (c/cmax)^1/3
    double cmax = M_PI/6.0;//for bcc cmax = 3^0.5*pi/8
    double R = std::pow(XX*XX*XX*cmax*h1*h1*h1/(4.0/3.0*M_PI),1.0/3.0);//x = (c/cmax)^1/3 c = Vs/V
    std::cout<<"R = "<<R<<std::endl;
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(0.0,0.0,2e-5);
    std::cout<<"nu "<<nu<<std::endl;
    std::cout<<"Nx "<<nx<<std::endl;
    std::cout<<"Ny "<<ny<<std::endl;
    std::cout<<"Nz "<<nz<<std::endl;
    LBM::Domain dom(D3Q19,MRT,nu,iVec3_t(nx,ny,nz),dx,dt);
    dom.Nproc = Nproc;
    dom.Sc = 0.0;
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.R = R;
    my_dat.u = 0.04;
    my_dat.g = 2e-5;
    my_dat.Tf = Tf;
    my_dat.K = 0.0;
    my_dat.nu = nu;
    my_dat.bbtype = bbtype;

    if(tau<0.53)
    {
        dom.S = 0, 1.19, 1.4, 0, 1.2, 0, 1.2, 0, 1.2, 1/tau, 1.4, 1/tau, 1.4, 1/tau, 1/tau, 1/tau, 1.98, 1.98, 1.98;
        dom.we = 0.0;
        dom.wej = -475.0/63.0;
        dom.wxx = 0;
    }

    //initial
    Initial(dom,rho0,v0,g0);
    std::fstream ifile("Nodes_sc.txt",std::ios::in);
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
            if(std::fabs(dot(post - pos, post - pos)-R*R)>1e-5) 
            {
                std::cout<<std::fabs(dot(post - pos, post - pos)-R*R)<<std::endl;
            }
            dom.points.push_back(post);
            dom.dS.push_back(ds);
		}
	}
    std::cout<<dom.dS[0]<<std::endl; 

    Tf = 2;
    double dtout = 1;
    char const * TheFileKey = "test_ibm";
    //solving
    dom.StartTime = std::clock();
    dom.StartSolve();
    double tout = 0;
    while(dom.Time<Tf)
    {
        
        if (dom.Time>=tout)
        {
            
            // String fn;
            // fn.Printf("%s_%d_%g_%04d", TheFileKey,bbtype, tau,dom.idx_out);
            // fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            // dom.WriteXDMF(fn.CStr());
            // dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        // (dom.*dom.ptr2collide)();
        
        
        dom.ApplyIBM3D(pos,R);
        dom.CollideMRTIBM();
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