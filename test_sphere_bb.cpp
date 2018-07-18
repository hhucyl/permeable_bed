#include "./lbm/Domain.h"
#include <fstream>

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
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"U"<<Util::_8s<<"K\n";
    }else{
        size_t index = 0;
        // double P1 = 0;
        // double num = 0;
        // for(size_t ix=0;ix<nx;++ix)
        // for(size_t iy=0;iy<ny;++iy)
        // {
        //     // if(fluid_dom.IsSolid[ix][iy][index]) continue;
        //     // num += 1;
        //     P1 += fluid_dom.Rho[ix][iy][index]/3.0;
        // }
        // P1 /= nx*ny;
        // index = nz -1;

        // double P2 = 0;
        // num = 0;
        // for(size_t ix=0;ix<nx;++ix)
        // for(size_t iy=0;iy<ny;++iy)
        // {
        //     // if(fluid_dom.IsSolid[ix][iy][index]) continue;
        //     // num +=1.0;
        //     P2 += fluid_dom.Rho[ix][iy][index]/3.0;
        // }
        // P2 /= nx*ny;
        double U = 0;
        double num = 0;
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
        // U += 0.5*dat.g;
        // std::cout<<dat.g<<" "<<dat.R<<std::endl;
        // double K = -(U*dat.nu)/((P2-P1)/((double)nz) - dat.g);
        double K = (U*dat.nu)/(dat.g);
        // double K = 8.0*dat.nu*nz*nz*U/(9*nz*(P1-P2));
        K = (6.0*M_PI*dat.R*K)/(nx*ny*nz);
        dat.oss_ss<<Util::_10_6<<fluid_dom.Time<<Util::_8s<<U<<Util::_8s<<K<<std::endl;
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
    size_t collidetype = 1;
    CollideMethod methodc = MRT;
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
    size_t Nproc = 8;
    size_t h = 100;
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double dx = 1.0;
    double dt = 1.0;
    double R = 10.0;
    Vec3_t pos(0.0,0.0,0.0);
    double nu = (0.6-0.5)/3.0;
    //nu = 1.0/30.0;
    LBM::Domain dom(D3Q19,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = 2e-5;
    my_dat.R = R;
    Vec3_t g0(0.0,0.0,my_dat.g);
    dom.Nproc = Nproc;


    std::fstream ifile("sphere.txt",std::ios::in);
	int N = 0;
	if(!ifile.fail())
	{
		ifile>>N;
		for(size_t i=0;i<N;++i)
		{
			double x;
			double y;
			double z;
    		ifile>>x>>y>>z;
            pos = x,y,z;
            dom.AddSphereQ(pos,R);
			
		}
	}          

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    //bounndary
    // dom.AddDisk(pos,R);
    //initial
    double rho = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0,g0);
    

    

    double Tf = 1e5;
    double dtout = 100;
    char const * TheFileKey = "test_sphere";
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
            // std::cout<<"--- Time = "<<dom.Time<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        (dom.*dom.ptr2collide)();
        dom.Stream();
        dom.BounceBackLIBB(false);
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    return 0;
}MECHSYS_CATCH