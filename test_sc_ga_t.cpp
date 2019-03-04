#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    int bbtype;
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
        fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
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

void addspheres1(LBM::Domain &dom, Vec3_t &pos,double R, double ddx)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    dom.AddSphereG(pos,R);
    pos = nx-1-ddx,0+ddx,0+ddx;
    dom.AddSphereG(pos,R);
    pos = 0+ddx,ny-1-ddx,0+ddx;
    dom.AddSphereG(pos,R);
    pos = nx-1-ddx, ny-1-ddx, 0+ddx;
    dom.AddSphereG(pos,R);
    pos = 0.0+ddx,0.0+ddx,nz-1-ddx;
    dom.AddSphereG(pos,R);
    pos = nx-1-ddx,0+ddx,nz-1-ddx;
    dom.AddSphereG(pos,R);
    pos = 0+ddx,ny-1-ddx,nz-1-ddx;
    dom.AddSphereG(pos,R);
    pos = nx-1-ddx, ny-1-ddx, nz-1-ddx;
    dom.AddSphereG(pos,R);  
}

void addspheres2(LBM::Domain &dom, Vec3_t &pos,double R, double ddx)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    dom.AddSphereG(pos,R);
}

//something to make compare simple
  
int main (int argc, char **argv) try
{
    size_t collidetype = 1;
    
    
    size_t Nproc = 8;
    size_t h = 30;
    int bbtype = -4;
    double tau = 0.6;
    if(argc>=2) bbtype = atoi(argv[1]); 
    if(argc>=3) h = atoi(argv[2]);
    if(argc>=4) tau = atof(argv[3]);     
    if(argc>=5) Nproc = atoi(argv[4]);
    if(argc>=6) collidetype = atoi(argv[5]);

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

    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double dx = 1.0;
    double dt = 1.0;
    double XX = 0.85;//x = (c/cmax)^1/3
    double cmax = M_PI/6.0;//for bcc cmax = 3^0.5*pi/8
    double R = std::pow(XX*XX*XX*cmax*nx*nx*nx/(4.0/3.0*M_PI),1.0/3.0);//x = (c/cmax)^1/3 c = Vs/V
    std::cout<<"R = "<<R<<std::endl;
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
    my_dat.g = 2e-5;
    my_dat.R = R;
    my_dat.K = 0;
    Vec3_t g0(0.0,0.0,my_dat.g);
    dom.Nproc = Nproc;
    if(tau<0.53)
    {
        dom.S = 0, 1.19, 1.4, 0, 1.2, 0, 1.2, 0, 1.2, 1/tau, 1.4, 1/tau, 1.4, 1/tau, 1/tau, 1/tau, 1.98, 1.98, 1.98;
        dom.we = 0.0;
        dom.wej = -475.0/63.0;
        dom.wxx = 0;
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
    
    addspheres2(dom,pos,R,ddx);
    

    double Tf = 1e6;
    my_dat.Tf = Tf;
    
    double dtout = 1e2;
    char const * TheFileKey = "test_sc";
    //solving
    dom.StartSolve();
    double tout = 0;
    while(dom.Time<Tf)
    {
        if (dom.Time>=tout)
        {
            
            // String fn;
            // fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            // dom.WriteXDMF(fn.CStr());
            // dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" "<<Tf<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        // dom.BoundaryGamma();
        // addspheres2(dom,pos,R,ddx);
        
        dom.CollideMRTPSM();
        // (dom.*dom.ptr2collide)();
        
        dom.Stream();
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    return 0;
}MECHSYS_CATCH