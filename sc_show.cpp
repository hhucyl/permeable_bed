//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

int main(int argc, char **argv) try
{
    size_t Nproc = 12; 
    if (argc>=2) Nproc=atoi(argv[1]);
    size_t h = 130;
    size_t nx = h;
    size_t ny = h;
    size_t nz = h;
    double nu = (1.6-0.5)/3.0;
    double dx = 1.0;
    double dt = 1.0;
    double Tf = 100000.0;
    double g  = 2e-8;
    Vec3_t g0(0.0,g,0.0);
    double XX = 0.85;//x = (c/cmax)^1/3
    double cmax = M_PI/6.0;//for bcc cmax = 3^0.5*pi/8
    double R = std::pow(XX*XX*XX*cmax*nx*nx*nx/(4.0/3.0*M_PI),1.0/3.0);//x = (c/cmax)^1/3 c = Vs/V
    std::cout<<"R = "<<R<<std::endl;
    if (argc>=3) nu =atof(argv[2]);
    

    LBM::Domain Dom(D3Q19, nu, iVec3_t(nx,ny,nz), dx, dt);
    Dom.Step     = 1;
    Dom.Sc       = 0.0;
    Dom.Alpha    = dx;
    

    
    Vec3_t pos(nx/2.0-1,ny/2.0-1,nz/2.0-1);
	
		
	
	
    Dom.AddSphere(-1, pos, R, 2.0);
    Dom.Particles[0]->FixVeloc();
    // Dom.Particles[i]->vzf = false;   
	

    
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
        Dom.Lat[0].Cells[i]->BForcef = g0;
    }

    //Solving
    Dom.Solve(100000,1000,NULL,NULL,"sc",true,Nproc);
}
MECHSYS_CATCH