#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double R;
    double Tff;
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
        fs.Printf("%s_%d_%g.out","Cd_LIBB_MRT",nx,dom.Tau);
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"Re"<<Util::_8s<<"Cd_ref1"<<Util::_8s<<"Cd_ref2"<<Util::_8s<<"Cd\n";
    }else{
        double U = 0;
        double num = 0;
        double Rho = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        for(size_t iz=0;iz<nz;++iz)
        {
            if(dom.IsSolid[ix][iy][iz]) continue;
            num += 1.0;
            Vec3_t e(0,1,0);
            U += dot(dom.Vel[ix][iy][iz],e);
            Rho += dom.Rho[ix][iy][iz];
        }
        U /= num;
        Rho /= num;
        double F = 0;
        num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        for(size_t iz=0;iz<nz;++iz)
        {
            F += dom.Flbm[ix][iy][iz](1);
        }
        double nu = (dom.Tau-0.5)/3.0;
        double Re = 2.0*dat.R*U/nu;
        double Cd = 2.0*F/(Rho*U*U*3.1415926*dat.R*dat.R);
        double Cd1 = 24.0/Re + 6.0/(1.0+std::sqrt(Re))+0.4;
        double Cd2 = 24.0/(9.06*9.06)*(9.06/std::sqrt(Re)+1)*(9.06/std::sqrt(Re)+1);
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<Re<<Util::_8s<<Cd1<<Util::_8s<<Cd2<<Util::_8s<<Cd<<std::endl;
        if(U >= 0.3) dom.Time = dat.Tff;
        
    }
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
    size_t bbtype = 3;    
    
    
    size_t Nproc = 12;
    size_t h = 100;
    double tau = 0.8;
    if(argc>=2) h = atoi(argv[1]); 
    if(argc>=3) tau = atof(argv[2]);
    if(argc>=4) Nproc = atoi(argv[3]); 
     
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double dx = 1.0;
    double dt = 1.0;
    double R = 5;
    double nu = (tau-0.5)/3.0;
    Vec3_t pos(49.5,48.5,49.5);
    // Vec3_t pos(200.5,299.5,0);
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(0.0,1e-6,0.0);
    
    LBM::Domain dom(D3Q19,methodc,nu,iVec3_t(nx,ny,nz),dx,dt);
    if(bbtype == 0)
    {
        dom.MethodB = SBB;
    }else if(bbtype == 1){
        dom.MethodB = LIBB;        
    }else if(bbtype == 2){
        dom.MethodB = QIBB;
    }else if(bbtype == 3){
        dom.MethodB = CLI;        
    }else{
        throw new Fatal("Collide Type is NOT RIGHT!!!!!");    
    }
    dom.SetBounceBack();

    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.R = R;
    dom.Nproc = Nproc;
    //debug mod
    // dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    //initial
    Initial(dom,rho0,v0,g0);
    dom.AddSphereQ(pos,R);

    double Tf = 1e2;
    my_dat.Tff = Tf;
    double dtout = 1;
    char const * TheFileKey = "test_carveb";
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
        (dom.*dom.ptr2bb)(true);
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    
    return 0;
}MECHSYS_CATCH