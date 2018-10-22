#include "./lbm/Domain.h"
#include <fstream>

struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double K;
    double Tf;
    size_t bbtype;
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
        //fs.Printf("%s.out","permeability");
        fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,fluid_dom.Tau);
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"U"<<Util::_8s<<"r"<<Util::_8s<<"K\n";
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
            Vec3_t e(0,1,0);
            U += dot(fluid_dom.Vel[ix][iy][0],e);
        }
        U /= nx*ny*nz;
        // U /= num;
        // U += 0.5*dat.g;
        // std::cout<<dat.g<<" "<<dat.R<<std::endl;
        // double K = -(U*dat.nu)/((P2-P1)/((double)nz) - dat.g);
        double K = (U*dat.nu)/(dat.g);
        // double K = 8.0*dat.nu*nz*nz*U/(9*nz*(P1-P2));
        // K = (6.0*M_PI*dat.R*K)/(nx*ny*nz);
		K /= nx*ny;
        double r = std::fabs(K-dat.K)/dat.K;
        dat.oss_ss<<Util::_10_6<<fluid_dom.Time<<Util::_8s<<U<<Util::_8s<<r<<Util::_8s<<K<<std::endl;
        if(r<1e-5&&U<0&&std::fabs(U)>1.0)
        {
            fluid_dom.Time = dat.Tf;
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

int main (int argc, char **argv) try
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
    int h = 200;
    size_t bbtype = 0;
    double tau = 2.0;
    if(argc>=2) bbtype = atoi(argv[1]); 
    if(argc>=3) h = atoi(argv[2]);
    if(argc>=4) tau = atof(argv[3]);     
    if(argc>=5) Nproc = atoi(argv[4]); 
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double dx = 1.0;
    double dt = 1.0;
    double R = 10.0;
    Vec3_t pos(0.0,0.0,0.0);
    double nu = (tau-0.5)/3.0;
    //nu = 1.0/30.0;
    LBM::Domain dom(D3Q19,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    if(bbtype == 0)
    {
        dom.MethodB = SBB;
    }else if(bbtype == 1){
        dom.MethodB = LIBB;        
    }else if(bbtype == 2){
        dom.MethodB = QIBB;
    }else if(bbtype == 3){
        dom.MethodB = MR;
    }else if(bbtype == 4){
        dom.MethodB = CLI;            
    }else{
        throw new Fatal("Collide Type is NOT RIGHT!!!!!");    
    }
    dom.SetBounceBack();
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = 2e-8;
    my_dat.R = R;
    my_dat.K= 0.0;
    my_dat.bbtype = bbtype;
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
        std::cout<<"size "<<h<<" R "<<R<<" Num "<<N<<std::endl;
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
    double dtout = 1e2;
    char const * TheFileKey = "test_sphere1";
    //solving
    dom.StartSolve();
    
    double tout = 0;
    while(dom.Time<Tf)
    {
        if (dom.Time>=tout)
        {
            
            //String fn;
            //fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            //dom.WriteXDMF(fn.CStr());
            //dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        (dom.*dom.ptr2collide)();
        // dom.CollideMRTMR();
        dom.Stream();
        // dom.BounceBackMR(false);
        (dom.*dom.ptr2bb)(false);
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    return 0;
}MECHSYS_CATCH
