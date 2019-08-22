#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    int bbtype;
    int collidetype;
    double K;
    double Tf;
};

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
        fs.Printf("%s_%d_%d_%d_%g.out","Permeability",dat.bbtype,dat.collidetype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"U"<<Util::_8s<<"r"<<Util::_8s<<"K\n";
    }else{
        size_t index = 0;
        // double P1 = 0;
        // double num = 0;
        // for(size_t ix=0;ix<nx;++ix)
        // for(size_t iy=0;iy<ny;++iy)
        // {
        //     // if(dom.IsSolid[ix][iy][index]) continue;
        //     // num += 1;
        //     P1 += dom.Rho[ix][iy][index]/3.0;
        // }
        // P1 /= nx*ny;
        // index = nz -1;

        // double P2 = 0;
        // num = 0;
        // for(size_t ix=0;ix<nx;++ix)
        // for(size_t iy=0;iy<ny;++iy)
        // {
        //     // if(dom.IsSolid[ix][iy][index]) continue;
        //     // num +=1.0;
        //     P2 += dom.Rho[ix][iy][index]/3.0;
        // }
        // P2 /= nx*ny;
        double U = 0;
        double num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        for(size_t iz=0;iz<nz;++iz)
        {
            // if(dom.IsSolid[ix][iy][iz]) continue;
            // num += 1.0;
            Vec3_t e(0,1,0);
            U += dot(dom.Vel[ix][iy][iz],e);
        }
        U /= nx*ny*nz;
        // U /= num;
        // U += 0.5*dat.g;
        // std::cout<<dat.g<<" "<<dat.R<<std::endl;
        // double K = -(U*dat.nu)/((P2-P1)/((double)nz) - dat.g);
        double K = (U*dat.nu)/(dat.g);
        // double K = 8.0*dat.nu*nz*nz*U/(9*nz*(P1-P2));
        K = K/(nx*ny);
        // K = (6.0*M_PI*dat.R*K)/(nx*ny*nz);
        double r = std::fabs(K-dat.K)/dat.K;
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<U<<Util::_8s<<r<<Util::_8s<<K<<std::endl;
        if(r<1e-5||U<0||std::fabs(U)>1.0)
        {
            dom.Time = dat.Tf;
        }
        dat.K = K;
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


void addspheres2(LBM::Domain &dom, Vec3_t &pos,double R, double *gammat, int nn)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    dom.AddSphereG(pos,R);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        gammat[ix+iy*nx+iz*nx*ny+nn*nx*ny*nz] = dom.Gamma[ix][iy][iz];
    }
}

//something to make compare simple
  
int main (int argc, char **argv) try
{
    int collidetype = 1;
    CollideMethod methodc = MRT;
    
    
    size_t Nproc = 8;
    size_t h = 50;
    int bbtype = -1;
    double tau = 0.6;
    if(argc>=2) bbtype = atoi(argv[1]); 
    if(argc>=3) h = atoi(argv[2]);
    if(argc>=4) tau = atof(argv[3]);     
    if(argc>=5) Nproc = atoi(argv[4]); 
    if(argc>=6) collidetype = atoi(argv[5]);
    if(collidetype == 0)
    {
        methodc = SRT;
        // ptr2collide = &LBM::Domain::CollideSRT;
    }else if(collidetype == 1){
        methodc = MRT;
        // ptr2collide = &LBM::Domain::CollideMRT;
        // ptr2meq = &LBM::Domain::MeqD2Q9;         
    }else{
        throw new Fatal("Collide Type is NOT RIGHT!!!!!");    
    }
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double dx = 1.0;
    double dt = 1.0;
    double R=0;
    double ddx = 0.0;
    // Vec3_t pos(ddx,ddx,ddx);
    Vec3_t pos(nx/2.0-1,ny/2.0-1,nz/2.0-1);
    // Vec3_t pos(200.5,299.5,0);
    double nu = (tau-0.5)/3.0;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D3Q19,methodc, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.bbtype = bbtype;
    my_dat.g = 2e-8;
    my_dat.R = R;
    my_dat.K = 0;
    my_dat.collidetype = collidetype;
    Vec3_t g0(0.0,my_dat.g,0.0);
    dom.Nproc = Nproc;

    char const *infilename = "spheres";
    String fn1;
    fn1.Printf("%s%d",infilename,h);
    fn1.append(".txt");
    std::fstream ifile(fn1.CStr(),std::ios::in);
	int N = 0;
    
	if(!ifile.fail())
	{
		ifile>>h>>R>>N;
        double *gammat = new double[nx*ny*nz*N];
        memset(gammat,0,nx*ny*nz*N); 
        std::cout<<"size "<<h<<" R "<<R<<" Num "<<N<<std::endl;
		for(size_t i=0;i<N;++i)
		{
			double x;
			double y;
			double z;
    		ifile>>x>>y>>z;
            pos = x,y,z;
            addspheres2(dom,pos,R,gammat,i);
			
		}
        
        for (size_t ix=0;ix<nx;ix++)
        for (size_t iy=0;iy<ny;iy++)
        for (size_t iz=0;iz<nz;iz++)
        {
            dom.Gamma[ix][iy][iz] = 0.0;
            for(size_t i=0; i<N; ++i)
            {
               dom.Gamma[ix][iy][iz] += gammat[ix+iy*nx+iz*nx*ny+i*nx*ny*nz];
            }
            dom.Gamma[ix][iy][iz] = std::min(1.0,dom.Gamma[ix][iy][iz]);
            
        }
        delete []gammat;
	}
    std::cout<<"Read Finish"<<std::endl;   

              

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    //bounndary
    // dom.AddDisk(pos,R);
    //initial
    double rho = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0,g0);
    

    

    double Tf = 1e6;
    my_dat.Tf = Tf;
    
    double dtout = 1e3;
    char const * TheFileKey = "test_spheres";
    //solving
    dom.StartSolve();
    double tout = 0;
    while(dom.Time<Tf)
    {
        if (dom.Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            dom.WriteXDMF(fn.CStr());
            dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" "<<Tf<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        dom.BoundaryGamma();
        // addspheres2(dom,pos,R,ddx);
        
        // dom.CollideSRT(ptr2meq);
        (dom.*dom.ptr2collide)();
        
        dom.Stream();
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    return 0;
}MECHSYS_CATCH
