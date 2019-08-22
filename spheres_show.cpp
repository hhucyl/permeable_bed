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
    if (argc>=3) nu =atof(argv[2]);
    

    LBM::Domain Dom(D3Q19, nu, iVec3_t(nx,ny,nz), dx, dt);
    Dom.Step     = 1;
    Dom.Sc       = 0.0;
    Dom.Alpha    = dx;
    // UserData dat;
    // Dom.UserData = &dat;
    // dat.Tf       = Tf;
    // dat.rhobc    = 1.0;
    // dat.g        = g;
    // dat.dx       = dx;
    // dat.Cs       = dx/(dt*sqrt(3.0));

    char const *infilename = "spheres";
    String fn1;
    fn1.Printf("%s%d",infilename,h);
    fn1.append(".txt");
    std::fstream ifile(fn1.CStr(),std::ios::in);
	int N = 0;
    double R=0;
    Vec3_t pos(0,0,0);
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
            Dom.AddSphere(-1, pos, R, 2.0);
            Dom.Particles[i]->FixVeloc();
            // Dom.Particles[i]->vzf = false;   
			
		}
        
	}
    std::cout<<"Read Finish"<<std::endl; 
    
    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
        Dom.Lat[0].Cells[i]->BForcef = g0;
    }

    //Solving
    Dom.Solve(5000,1000,NULL,NULL,"ss",true,Nproc);
}
MECHSYS_CATCH