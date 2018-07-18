#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double qt;
    double Tff;
    double uw;
    double g;
    size_t collidetype;
    size_t bbtype;
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
}

void Setup(LBM::Domain &dom,  void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
    for(size_t ix=0; ix<nx; ix++)
    for(size_t iz=0; iz<nz; iz++)
    {
        dom.VelP[ix][ny-1][iz] = dat.uw,0.0,0.0;
        // dom.VelP[ix][ny-2][iz] = dat.uw,0.0,0.0;
    }
   
}

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    double dx = dom.dx;
    double nu = dom.Nu;
    if(dom.Time <1e-6)
    {
        String fs;
        // fs.Printf("%s_%d_%g.out","R",dat.bbtype,dom.Tau);
        fs.Printf("%s.out","R");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"u_ref"<<Util::_8s<<"u"<<Util::_8s<<"r\n";
    }else{
        double u;
        double ur;
        double y0 = 0.5;
        double y1 = ny-1+dat.qt;
        double H = (y1-y0);
        double r = 0;
        size_t indexx = nx/2;
        size_t indexz = nz/2;
        for(size_t iy=0;iy<ny-1;++iy)
        {
            if(dom.IsSolid[indexx][iy][indexz]) continue;
            u = dom.Vel[indexx][iy][indexz](0);
            double yr = (double)iy-0.5;
            ur = -dat.g/(2.0*nu)*yr*yr + (dat.g/(2.0*nu)*H+dat.uw/H)*yr;
            // std::cout<<H<<" "<<iy<<" "<<yr<<" "<<u<<" "<<ur<<std::endl;
            r += fabs((u-ur)/ur);
        }
        r /= ny-1;
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<ur<<Util::_8s<<u<<Util::_8s<<r<<std::endl;
        // if(U >= 0.3) dom.Time = dat.Tff;
        
    }
}

using namespace std;

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
    size_t bbtype = 3;    
    size_t Nproc = 10;
    size_t h = 50;
    double tau = 0.8;
    // if(argc>=1) bbtype = atoi(argv[1]); 
    // if(argc>=2) tau = atof(argv[2]);
    // if(argc>=3) Nproc = atoi(argv[3]); 

    
    double qt = 0.25; 
    size_t nx = 2.0*h;
    size_t ny = h;
    if(qt<0.5)
    {
        ny += 0;
    }else{
        ny += 1;
    }
    size_t nz = 10;
    double dx = 1.0;
    double dt = 1.0;
    double nu = (tau-0.5)/3.0;
    double uw = 0.01;
    // Vec3_t pos(200.5,299.5,0);
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(1.0e-6,0.0,0.0);
    
    LBM::Domain dom(D3Q19,methodc,nu,iVec3_t(nx,ny,nz),dx,dt);
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
    my_dat.qt = qt;
    my_dat.g = g0(0);
    my_dat.collidetype = collidetype;
    my_dat.bbtype = bbtype;
    my_dat.uw = uw;
    dom.Nproc = Nproc;
    //debug mod
    // dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    //boundary
   
    for(size_t ix=0; ix<nx; ix++)
    for(size_t iz=0; iz<nz; iz++)
    {
        dom.IsSolid[ix][0][iz] = true;
        // dom.Gamma[ix][ny-1][iz] = 1.0;
        if(qt<0.5)
        {
            dom.Gamma[ix][ny-1][iz] = 0.5 - qt;        
        }else{
            dom.Gamma[ix][ny-1][iz] = 1.5 - qt;
            
        }
        // dom.VelP[ix][ny-2][iz] = uw,0.0,0.0;
        dom.VelP[ix][ny-1][0] = uw,0.0,0.0;
    }
   
    //initial
    Initial(dom,rho0,v0,g0);


    double Tf = 1e5;
    my_dat.Tff = Tf;
    double dtout = 1e2;
    char const * TheFileKey = "test_line";
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
            fn.Printf("%s_%04d", TheFileKey,dom.idx_out);
            
            dom.WriteXDMF(fn.CStr());
            dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        // (dom.*dom.ptr2collide)();
        Setup(dom,&my_dat);
        dom.BoundaryGamma();
        dom.CollideMRT();
        dom.Stream();
        // dom.BounceBack(false);
        // dom.BounceBackMR(false);
        // (dom.*dom.ptr2bb)(true);
        dom.CalcProps();
        dom.Time += 1;
    }
    // dom.EndSolve();
    
    dom.EndTime = std::clock();
    double ttime = (double)(dom.EndTime - dom.StartTime)/CLOCKS_PER_SEC;
    printf("\033[01;34m---Elapsed Time = %f s---\033[0m\n",ttime);
    
    return 0;
}MECHSYS_CATCH