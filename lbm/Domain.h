

#ifndef MECHSYS_LBM_DOMAIN1_H
#define MECHSYS_LBM_DOMAIN1_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#ifdef USE_OCL
#include <mechsys/oclaux/cl.hpp>
#endif

// Std lib
#ifdef USE_OMP
#include <omp.h>
#endif

//STD
#include<iostream>

// Mechsys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/numstreams.h>

enum LBMethod
{
    D2Q5,     ///< 2D 5 velocities
    D2Q9,     ///< 2D 9 velocities
    D3Q15,    ///< 3D 15 velocities
    D3Q19,    ///< 3D 19 velocities
    //D3Q27     ///< 3D 27 velocities
};

enum CollideMethod
{
    SRT,
    MRT
};

enum BounceBackMethod
{
    SBB,
    LIBB,
    QIBB
};



namespace LBM
{



inline size_t Pt2idx(iVec3_t iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}

class Domain
{
public:
	static const double   WEIGHTSD2Q5   [ 5]; ///< Weights for the equilibrium distribution functions (D2Q5)
	static const double   WEIGHTSD2Q9   [ 9]; ///< Weights for the equilibrium distribution functions (D2Q9)
	static const double   WEIGHTSD3Q15  [15]; ///< Weights for the equilibrium distribution functions (D3Q15)
	static const double   WEIGHTSD3Q19  [19]; ///< Weights for the equilibrium distribution functions (D3Q19)
	//static const double   WEIGHTSD3Q27  [27]; ///< Weights for the equilibrium distribution functions (D3Q27)
	static const Vec3_t   LVELOCD2Q5    [ 5]; ///< Local velocities (D2Q5) 
	static const Vec3_t   LVELOCD2Q9    [ 9]; ///< Local velocities (D2Q9) 
	static const Vec3_t   LVELOCD3Q15   [15]; ///< Local velocities (D3Q15)
	static const Vec3_t   LVELOCD3Q19   [19]; ///< Local velocities (D3Q19)
	//static const Vec3_t   LVELOCD3Q27   [27]; ///< Local velocities (D3Q27)
	static const size_t   OPPOSITED2Q5  [ 5]; ///< Opposite directions (D2Q5) 
	static const size_t   OPPOSITED2Q9  [ 9]; ///< Opposite directions (D2Q9) 
	static const size_t   OPPOSITED3Q15 [15]; ///< Opposite directions (D3Q15)
	static const size_t   OPPOSITED3Q19 [19]; ///< Opposite directions (D3Q19)
	//static const size_t   OPPOSITED3Q27 [27]; ///< Opposite directions (D3Q27)
    static const double   MD2Q5       [5][5]; ///< MRT transformation matrix (D2Q5)
    static const double   MD2Q9       [9][9]; ///< MRT transformation matrix (D2Q9)
    static const double   MD3Q15    [15][15]; ///< MRT transformation matrix (D3Q15)
    static const double   MID3Q15    [15][15]; ///< MRT transformation matrix (D3Q15)
    static const double   MD3Q19    [19][19]; ///< MRT transformation matrix (D3Q19)
    //static const size_t   MD3Q27    [27][27]; ///< MRT transformation matrix (D3Q27)
    //static const double   SD2Q5          [5]; ///< MRT relaxation time vector (D2Q5)
    //static const double   SD2Q9          [9]; ///< MRT relaxation time vector (D2Q9)
    //static const double   SD3Q15        [15]; ///< MRT relaxation time vector (D3Q15)
    //static const double   SD3Q19        [19]; ///< MRT relaxation time vector (D3Q19)
    //static       double   SD3Q19        [27]; ///< MRT relaxation time vector (D3Q27)
    
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);
    typedef void (LBM::Domain::*ptFM_t)  (double *m, double rho, Vec3_t &vel); 
    typedef void (LBM::Domain::*ptFC_t)  (ptFM_t ptr);
    typedef void (LBM::Domain::*ptFB_t)  ();
    
   
    
    //Constructors
    

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    CollideMethod         MethodC, //< Type of collision
    BounceBackMethod      MethodB, //< Type of bounceback
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Methods
    
    void   CollideMRT(ptFM_t ptr);///< The collide step of LBM with MRT                                  
    void   MeqD2Q9(double *m, double rho, Vec3_t &vel);                      
    void   MeqD3Q15(double *m, double rho, Vec3_t &vel);             
    void   MeqD3Q19(double *m, double rho, Vec3_t &vel);             
    
                           ///< The collide step of LBM with MRT
    void   CollideSRT(ptFM_t ptr);                                                           ///< The collide step of LBM for single component simulations
    void   Stream();
    void   BounceBack();
    void   BounceBackLIBB();
    void   CalcProps();
    void   Initialize(iVec3_t idx, double Rho, Vec3_t & Vel);           ///< Initialize each cell with a given density and velocity
    double Feq(size_t k, double Rho, Vec3_t & Vel);                               ///< The equilibrium function
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics
    
    //Writing Methods
    void WriteXDMF         (char const * FileKey);                                ///< Write the domain data in xdmf file
    void AddDiskQ(Vec3_t &pos, double R);
    void CalcForceQ();
    
    
    
    #ifdef USE_OMP
    omp_lock_t      lck;                      ///< to protect variables in multithreading
    #endif

    //Data
    double ****q;
    double ****Omeis;
    double **** F;                           ///< The array containing the individual functions with the order of the lattice, the x,y,z coordinates and the order of the function.
    double **** Ftemp;                       ///< A similar array to hold provitional data
    bool   ***  IsSolid;                     ///< An array of bools with an identifier to see if the cell is a solid cell
    Vec3_t ***  Vel;                         ///< The fluid velocities
    Vec3_t ***  BForce;                      ///< Body Force for each cell
    double ***  Rho;                         ///< The fluid densities
    double      Tau;                         ///< The characteristic time of the lattice
                          ///< The mixing constant for multicomponent simulations
    size_t const * Op;                        ///< An array containing the indexes of the opposite direction for bounce back conditions
    double const *  W;                        ///< An array with the direction weights
    double *     EEk;                         ///< Diadic product of the velocity vectors
    Vec3_t const *  C;                        ///< The array of lattice velocities
    Mat_t        M;                           ///< Transformation matrix to momentum space for MRT calculations
    Mat_t        Minv;                        ///< Inverse Transformation matrix to momentum space for MRT calculations
    Vec_t        S;                           ///< Vector of relaxation times for MRT
    size_t       Nneigh;                      ///< Number of Neighbors, depends on the scheme
    double       dt;                          ///< Time Step
    double       dx;                          ///< Grid size
    double       Cs;                          ///< Lattice Velocity
    bool         IsFirstTime;                 ///< Bool variable checking if it is the first time function Setup is called
    iVec3_t      Ndim;                        ///< Lattice Dimensions
    size_t       Ncells;                      ///< Number of cells
    size_t       Nproc;                       ///< Number of processors for openmp
    size_t       idx_out;                     ///< The discrete time step for output
    String       FileKey;                     ///< File Key for output files
    void *       UserData;                    ///< User Data
    size_t       Step;                        ///< Lenght of averaging cube to save data
    double       Time;                        ///< Simulation time variable
    size_t       Nl;                          ///< Number of lattices (fluids)
    double       Sc;                          ///< Smagorinsky constant
    LBMethod     Method;
    CollideMethod MethodC;
    BounceBackMethod MethodB;
    bool IsF;
    bool IsFt;
    bool Isq;
    Vec3_t ***VelP;
    Vec3_t ***Flbm;
    double ***Gamma;

};



inline Domain::Domain(LBMethod TheMethod, CollideMethod TheMethodC, BounceBackMethod TheMethodB, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Array<double> nu(1);
    nu[0] = Thenu;
    
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (TheNdim(2) >1&&(TheMethod==D2Q9 ||TheMethod==D2Q5 ))  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (TheNdim(2)==1&&(TheMethod==D3Q15||TheMethod==D3Q19))  throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
    IsF = false;
    IsFt = false;
    Isq = false;
    Time        = 0.0;
    dt          = Thedt;
    dx          = Thedx;
    Cs          = dx/dt;
    Step        = 1;
    Sc          = 0.17;
    Nl          = 1;
    Ndim        = TheNdim;
    Ncells      = Ndim(0)*Ndim(1)*Ndim(2);
    IsFirstTime = true;
    Method = TheMethod;
    MethodC = TheMethodC;
    MethodB = TheMethodB;
    Tau         = 3.0*nu[0]*dt/(dx*dx)+0.5;
    if (TheMethod==D2Q5)
    {
        Nneigh = 5;
        W      = WEIGHTSD2Q5;
        C      = LVELOCD2Q5;
        Op     = OPPOSITED2Q5;
    }
    if (TheMethod==D2Q9)
    {
        Nneigh = 9;
        W      = WEIGHTSD2Q9;
        C      = LVELOCD2Q9;
        Op     = OPPOSITED2Q9;
        if(TheMethodC == MRT)
        {
            M.Resize(Nneigh,Nneigh);
            Minv.Resize(Nneigh,Nneigh);
            for (size_t n=0;n<Nneigh;n++)
            for (size_t m=0;m<Nneigh;m++)
            {
                M(n,m) = MD2Q9[n][m];
                // Minv(n,m) = MID3Q15[n][m];
            }
            Inv(M,Minv);
            
            double s = 1.0/Tau;
            S.Resize(Nneigh);
            S = 1.0,1.4,1.4,1.0,1.2,1.0,1.2,s,s;
        }
    }
    if (TheMethod==D3Q15)
    {
        Nneigh = 15;
        W      = WEIGHTSD3Q15;
        C      = LVELOCD3Q15;
        Op     = OPPOSITED3Q15;
        if(TheMethodC == MRT)
        {
            M.Resize(Nneigh,Nneigh);
            Minv.Resize(Nneigh,Nneigh);
            for (size_t n=0;n<Nneigh;n++)
            for (size_t m=0;m<Nneigh;m++)
            {
                M(n,m) = MD3Q15[n][m];
                // Minv(n,m) = MID3Q15[n][m];
            }
            Inv(M,Minv);
            double tau = 3.0*Thenu*dt/(dx*dx)+0.5;
            double s   = 8.0*(2.0-1.0/tau)/(8.0-1.0/tau);
            S.Resize(Nneigh);
            // S = 0.0,1.0/tau,1.0/tau,0.0,s,0.0,s,0.0,s,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,s;
            S = 0.0,1.6,1.2,0.0,1.6,0.0,1.6,0.0,1.6,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.2;
        }
    }
    if (TheMethod==D3Q19)
    {
        Nneigh = 19;
        W      = WEIGHTSD3Q19;
        C      = LVELOCD3Q19;
        Op     = OPPOSITED3Q19;
        if(TheMethodC == MRT)
        {
            M.Resize(Nneigh,Nneigh);
            Minv.Resize(Nneigh,Nneigh);
            for (size_t n=0;n<Nneigh;n++)
            for (size_t m=0;m<Nneigh;m++)
            {
                M(n,m) = MD3Q19[n][m];
                // Minv(n,m) = MID3Q15[n][m];
            }
            Inv(M,Minv);
            double tau = 3.0*Thenu*dt/(dx*dx)+0.5;
            double s   = 8.0*(2.0-1.0/tau)/(8.0-1.0/tau);
            S.Resize(Nneigh);
            // S = 0.0,1.0/tau,1.0/tau,0.0,s,0.0,s,0.0,s,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,s;
            S = 0.0,s,s,0.0,s,0.0,s,0.0,s,1.0/tau,s,1.0/tau,s,1.0/tau,1.0/tau,1.0/tau,s,s,s;
        }
    }

    
    q           = new double *** [Ndim(0)];
	Omeis       = new double *** [Ndim(0)];
    VelP        = new Vec3_t **  [Ndim(0)];
    Flbm        = new Vec3_t **  [Ndim(0)];
    F           = new double *** [Ndim(0)];
    Ftemp       = new double *** [Ndim(0)];
    Vel         = new Vec3_t **  [Ndim(0)];
    BForce      = new Vec3_t **  [Ndim(0)];
    Rho         = new double **  [Ndim(0)];
    IsSolid     = new bool   **  [Ndim(0)];
    Gamma       = new double **  [Ndim(0)];
    for (size_t nx=0;nx<Ndim(0);nx++)
    {
        F       [nx]    = new double ** [Ndim(1)];
        Ftemp   [nx]    = new double ** [Ndim(1)];
        Vel     [nx]    = new Vec3_t *  [Ndim(1)];
        BForce  [nx]    = new Vec3_t *  [Ndim(1)];
        Rho     [nx]    = new double *  [Ndim(1)];
        IsSolid [nx]    = new bool   *  [Ndim(1)];
        q       [nx]    = new double ** [Ndim(1)];
	    Omeis   [nx]    = new double ** [Ndim(1)];
        VelP    [nx]    = new Vec3_t *  [Ndim(1)];
        Flbm    [nx]    = new Vec3_t *  [Ndim(1)];
        Gamma   [nx]    = new double *  [Ndim(1)];
        for (size_t ny=0;ny<Ndim(1);ny++)
        {
            F       [nx][ny]    = new double * [Ndim(2)];
            Ftemp   [nx][ny]    = new double * [Ndim(2)];
            Vel     [nx][ny]    = new Vec3_t   [Ndim(2)];
            BForce  [nx][ny]    = new Vec3_t   [Ndim(2)];
            Rho     [nx][ny]    = new double   [Ndim(2)];
            IsSolid [nx][ny]    = new bool     [Ndim(2)];
            q       [nx][ny]    = new double * [Ndim(2)];
	        Omeis   [nx][ny]    = new double * [Ndim(2)];
            VelP    [nx][ny]    = new Vec3_t   [Ndim(2)];
            Flbm    [nx][ny]    = new Vec3_t   [Ndim(2)];
            Gamma   [nx][ny]    = new double   [Ndim(2)];        
            for (size_t nz=0;nz<Ndim(2);nz++)
            {
                F    [nx][ny][nz]    = new double [Nneigh];
                Ftemp[nx][ny][nz]    = new double [Nneigh];
                q    [nx][ny][nz]    = new double [Nneigh];
	            Omeis[nx][ny][nz]    = new double [Nneigh];
                IsSolid[nx][ny][nz]  = false;
                VelP[nx][ny][nz] = 0.0, 0.0, 0.0;
                Flbm[nx][ny][nz] = 0.0, 0.0, 0.0;
                BForce[nx][ny][nz] = 0.0, 0.0, 0.0;
                Gamma[nx][ny][nz] = 0.0;
                for (size_t nn=0;nn<Nneigh;nn++)
                {
                    F    [nx][ny][nz][nn] = 0.0;
                    Ftemp[nx][ny][nz][nn] = 0.0;
                    q    [nx][ny][nz][nn] = -1.0;
                    Omeis[nx][ny][nz][nn] = 0.0;
                }
            }
        }
    }
    


    EEk = new double [Nneigh];
    for (size_t k=0;k<Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(C[k][n]*C[k][m]);
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Nl*Ncells,TERM_RST);
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

#include "Output.h"
#include "AddSolidBoundary.h"

inline void Domain::CalcForceQ()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx;ix++)
	for(size_t iy=0; iy<ny;iy++)
	for(size_t iz=0; iz<nz;iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            // size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            // size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            // size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k];
        }
        Flbm[ix][iy][iz] = Flbmt;
    }

}

inline double Domain::Feq(size_t k, double Rho, Vec3_t & V)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    return W[k]*Rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

inline void Domain::MeqD2Q9(double *m, double rho, Vec3_t &vel)
{
    m[0] = S(0)*(m[0] - rho);
    m[1] = S(1)*(m[1] - rho*(-2.0 + 3.0*dot(vel,vel)));
    m[2] = S(2)*(m[2] - rho*( 1.0 - 3.0*dot(vel,vel)));
    m[3] = S(3)*(m[3] - rho*vel(0));
    m[4] = S(4)*(m[4] + rho*vel(0));
    m[5] = S(5)*(m[5] - rho*vel(1));
    m[6] = S(6)*(m[6] + rho*vel(1));
    m[7] = S(7)*(m[7] - rho*(vel(0)*vel(0)-vel(1)*vel(1)));
    m[8] = S(8)*(m[8] - rho*(vel(0)*vel(1)));
}

inline void Domain::MeqD3Q15(double *m, double rho, Vec3_t &vel)
{
    m[ 0] = S(0)*(m[0] - rho); 
    m[ 1] = S( 1)*(m[ 1] + rho - rho*dot(vel,vel));
    m[ 2] = S( 2)*(m[ 2] + rho);
    m[ 3] = 0.0;
    m[ 4] = S( 4)*(m[ 4] + 7.0/3.0*rho*vel(0)); 
    m[ 5] = 0.0;
    m[ 6] = S( 6)*(m[ 6] + 7.0/3.0*rho*vel(1)); 
    m[ 7] = 0.0;
    m[ 8] = S( 8)*(m[ 8] + 7.0/3.0*rho*vel(2)); 
    m[ 9] = S( 9)*(m[ 9] - rho*(3.0*vel(0)*vel(0)-dot(vel,vel)));
    m[10] = S(10)*(m[10] - rho*(vel(1)*vel(1)-vel(2)*vel(2)));
    m[11] = S(11)*(m[11] - rho*(vel(0)*vel(1)));
    m[12] = S(12)*(m[12] - rho*(vel(1)*vel(2)));
    m[13] = S(13)*(m[13] - rho*(vel(0)*vel(2)));
    m[14] = S(14)* m[14];
}

inline void Domain::MeqD3Q19(double *m, double rho, Vec3_t &vel)
{
    m[ 0] = S( 0)*( m[ 0] - rho );
    m[ 1] = S( 1)*( m[ 1] - (-11.0*rho + 19.0 * dot(vel,vel)) );
    m[ 2] = S( 2)*( m[ 2] - (3.0*rho - 11.0/2.0*rho*dot(vel,vel)) );
    m[ 3] = 0.0;
    m[ 4] = S( 4)*( m[ 4] - (-2.0/3.0*rho*vel(0)) );
    m[ 5] = 0.0;
    m[ 6] = S( 6)*( m[ 6] - (-2.0/3.0*rho*vel(1)) );
    m[ 7] = 0.0;
    m[ 8] = S( 8)*( m[ 8] - (-2.0/3.0*rho*vel(2)) );
    m[ 9] = S( 9)*( m[ 9] - (2.0*rho*vel(0)*vel(0)-rho*vel(1)*vel(1)-rho*vel(2)*vel(2)) );
    m[10] = S(10)*( m[10] - (-0.5*(2.0*rho*vel(0)*vel(0)-rho*vel(1)*vel(1)-rho*vel(2)*vel(2))) );
    m[11] = S(11)*( m[11] - rho*(vel(1)*vel(1)-vel(2)*vel(2)) );
    m[12] = S(12)*( m[12] - (-0.5*rho*(vel(1)*vel(1)-vel(2)*vel(2))) );
    m[13] = S(13)*( m[13] - rho*vel(0)*vel(1) );
    m[14] = S(14)*( m[14] - rho*vel(1)*vel(2) );
    m[15] = S(15)*( m[15] - rho*vel(0)*vel(2) );
    m[16] = S(16)*  m[16];
    m[17] = S(17)*  m[17];
    m[18] = S(18)*  m[18];
}

inline void Domain::Initialize(iVec3_t idx, double TheRho, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);

    BForce[ix][iy][iz] = OrthoSys::O;

    for (size_t k=0;k<Nneigh;k++)
    {
        F[ix][iy][iz][k] = Feq(k,TheRho,TheVel);
    }

    if (!IsSolid[ix][iy][iz])
    {
        Vel[ix][iy][iz] = TheVel;
        Rho[ix][iy][iz] = TheRho;
    }
    else
    {
        Vel[ix][iy][iz] = OrthoSys::O;
        Rho[ix][iy][iz] = 0.0;
    }
}


inline void Domain::CollideSRT(ptFM_t ptr)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            double NonEq[Nneigh];
            // double Q = 0.0;
            double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz];//+dt*tau*BForce[ix][iy][iz]/rho;
            double VdotV = dot(vel,vel);
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,C[k]);
                double Feq   = W[k]*rho*(1.0 + 3.0*VdotC + 4.5*VdotC*VdotC - 1.5*VdotV);
                NonEq[k] = F[ix][iy][iz][k] - Feq;
            //     Q +=  NonEq[k]*NonEq[k]*EEk[k];
            }
            // Q = sqrt(2.0*Q);
            // tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));

            // bool valid = true;
            // double alpha = 1.0;
            // while (valid)
            // {
                // valid = false;
                for (size_t k=0;k<Nneigh;k++)
                {
                    double ForceTerm = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                    Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - NonEq[k]/tau + ForceTerm;
                    // if (Ftemp[ix][iy][iz][k]<-1.0e-12)
                    // {
                    //     //std::cout << Ftemp[ix][iy][iz][k] << std::endl;
                    //     double temp =  tau*F[ix][iy][iz][k]/NonEq[k];
                    //     if (temp<alpha) alpha = temp;
                    //     valid = true;
                    // }
                    // if (std::isnan(Ftemp[ix][iy][iz][k]))
                    // {
                    //     std::cout << "CollideSC: Nan found, resetting" << std::endl;
                    //     std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                    //     throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    // }
                }
            // }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }
}

inline void Domain::CollideMRT(ptFM_t ptr)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            //double NonEq[Nneigh];
            //double Q = 0.0;
            //double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz];
            // Vec3_t vel = Vel[ix][iy][iz]+0.5*dt*BForce[ix][iy][iz];
            //double VdotV = dot(vel,vel);
            //for (size_t k=0;k<Nneigh;k++)
            //{
                //double VdotC = dot(vel,C[k]);
                //double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                //NonEq[k] = F[ix][iy][iz][k] - Feq;
                //Q +=  NonEq[k]*NonEq[k]*EEk[k];
            //}
            //Q = sqrt(2.0*Q);
            //tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));
            double *f = F[ix][iy][iz];
            double *m = Ftemp[ix][iy][iz];
            double fneq[Nneigh];
            int n=Nneigh,mm=1;
            double a = 1.0,b = 0.0;
            dgemv_("N",&n,&n,&a,M.data,&n,f,&mm,&b,m,&mm);
            // for(size_t i=0; i<Nneigh; i++)
            // {
            //     m[i] = 0.0;
            //     for(size_t j=0; j<Nneigh; j++)
            //     {
            //         m[i] += M(i,j)*f[j];                   
            //     }
            // }
            
            
            (this->*ptr)(m,rho,vel);
            // meqD3Q15(m,rho,vel);
            
            dgemv_("N",&n,&n,&a,Minv.data,&n,m,&mm,&b,fneq,&mm);
            // for(size_t i=0; i<Nneigh; i++)
            // {
            //     fneq[i] = 0.0;
            //     for(size_t j=0; j<Nneigh; j++)
            //     {
            //         fneq[i] += Minv(i,j)*m[j];                   
            //     }
            // }
            for (size_t k=0; k<Nneigh; k++)
            {
                double ForceTerm = 1.0*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                // Vec3_t BFt(0.0, 0.0, 0.0);
                // BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
                // double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,BForce[ix][iy][iz]); 
                
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - fneq[k] + ForceTerm;
                // if(ix == 1 &&iy==10&&iz ==0)
                // {
                //     std::cout<<"F "<<k<<"="<<F[ix][iy][iz][k]<<std::endl;
                //     std::cout<<"Ft "<<k<<"="<<Ftemp[ix][iy][iz][k]<<std::endl;
                // }
            }

            //bool valid = true;
            //double alpha = 1.0, alphat = 1.0;
            //while (valid)
            //{
                //alpha = alphat;
                //valid = false;
              //  for (size_t k=0;k<Nneigh;k++)
            //    {
          //          Vec3_t BFt(0.0, 0.0, 0.0);
        //            BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
      //              double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,BForce[ix][iy][iz]); 
                    //Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - alpha*(fneq[k]-ForceTerm);
                    //if (Ftemp[ix][iy][iz][k]<-1.0e-12)
                    //{
                        //std::cout << Ftemp[ix][iy][iz][k] << std::endl;
                      //  valid = true;
                    //    double temp =  F[ix][iy][iz][k]/(fneq[k]-ForceTerm);
                  //      if (temp<alphat) alphat = temp;
                //    }
                    //if (std::isnan(Ftemp[ix][iy][iz][k]))
                    //{
                        //std::cout << "CollideSC: Nan found, resetting" << std::endl;
                        //std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                        //throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    //}
              //  }
            //}
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }

}

inline void Domain::Stream()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);
            F[nix][niy][niz][k] = Ftemp[ix][iy][iz][k] ;
        }
    }

}

inline void Domain::CalcProps()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        //BForce[ix][iy][iz] = OrthoSys::O;
        Vel   [ix][iy][iz] = OrthoSys::O;
        Rho   [ix][iy][iz] = 0.0;
        // if (!IsSolid[ix][iy][iz])
        // {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[ix][iy][iz] +=  F[ix][iy][iz][k];
                Vel[ix][iy][iz] +=  F[ix][iy][iz][k]*C[k];
            }
            Vel[ix][iy][iz] *= Cs/Rho[ix][iy][iz];
        // }
    }
}

inline void Domain::BounceBack()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        

        for(size_t k=0; k<Nneigh; k++)
        {
            size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
            
            if(IsSolid[oix][oiy][oiz])
            {   
                F[ix][iy][iz][k] = Ftemp[ix][iy][iz][Op[k]];
                
            }
            
        }

	    
    }
}

inline void Domain::BounceBackLIBB()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        double Temp[Nneigh];
        double Temp1[Nneigh];
        for(size_t k=0; k<Nneigh; k++)
        {
            Temp[k] = Ftemp[ix][iy][iz][k];
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)nx)%nx;
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)ny)%ny;
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)nz)%nz;
            Temp1[k] = Ftemp[nix][niy][niz][k];
        }
        
        
        for(size_t k=0; k<Nneigh; k++)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
            
            if(IsSolid[oix][oiy][oiz])
            {
                // Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k];
                if(qt<0.5)
                {
                    // if(IsSolid[nix][niy][niz]) qt = 0.5;
                    // F[ix][iy][iz][k] = 2.0*qt*Temp[Op[k]]+(1-2.0*qt)*F[ix][iy][iz][Op[k]];                 
                    F[ix][iy][iz][k] = 2.0*qt*Temp[Op[k]]+(1-2.0*qt)*Temp1[Op[k]];                 
                }else{
                    F[ix][iy][iz][k] = (1.0/(2.0*qt))*Temp[Op[k]]+(1-1.0/(2.0*qt))*Temp[k];                                                             
                }

            }
            
            
        }
        

        

    }   
}



inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{
    


    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;
    ptFM_t ptr2meq = NULL;
    ptFC_t ptr2collide = NULL;
    ptFB_t ptr2bounceback = NULL;
    ptFB_t ptr2calcforce = NULL;
    if(MethodC == MRT)
    {
        if(Method == D2Q9)
            ptr2meq = &LBM::Domain::MeqD2Q9;        
        if(Method == D3Q15)
            ptr2meq = &LBM::Domain::MeqD3Q15;
        if(Method == D3Q19)
            ptr2meq = &LBM::Domain::MeqD3Q19;

        ptr2collide = &LBM::Domain::CollideMRT;
    }else if(MethodC == SRT)
    {
        ptr2collide = &LBM::Domain::CollideSRT;
    }else{
        throw new Fatal("Collide Type is NOT RIGHT!!!!!");
    }
    if(MethodB == SBB)
    {
        ptr2bounceback = &LBM::Domain::BounceBack;
    }else if(MethodB == LIBB)
    {
        ptr2bounceback = &LBM::Domain::BounceBackLIBB; 
        ptr2calcforce = &LBM::Domain::CalcForceQ;       
    }

    if(ptr2collide == NULL)
    {
        throw new Fatal("Collsion is not initialized, check the pointer");
    }else{
        if(MethodC == MRT)
        {
            if(ptr2meq == NULL)
                throw new Fatal("Matrix for MRT is not initialized, check the pointer");
            
        }
    }
    if(ptr2bounceback == NULL)
    {
        printf("\n%s--- Warning pointer for bounceback is not assigned, now using 1st order SBB ----%s\n",TERM_CLR1    , TERM_RST);        

    }

    
    // if()
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving -------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    
    printf("%s  Tau of Lattice                  =  %g%s\n"       ,TERM_CLR2, Tau                             , TERM_RST);
    
  
    

    
    
    //std::cout << "2" << std::endl;

    double tout = Time;
    while (Time < Tf)
    {
        if (Time >= tout)
        {
            //std::cout << "3" << std::endl;
            #ifdef USE_OCL
            DnLoadDevice();
            #endif
            //std::cout << "4" << std::endl;
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    #else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
                std::cout<<"-------------Time = "<<Time<<" ----------------"<<std::endl;
            }
            tout += dtOut;
            idx_out++;
        }

        
        
        // CollideMRT(ptr2meq);
        (this->*ptr2collide)(ptr2meq);
        
        Stream();
        (this->*ptr2bounceback)();
        if (ptr2calcforce!=NULL)(this->*ptr2calcforce)();
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        
        CalcProps();
        
        

        Time += dt;
        //std::cout << Time << std::endl;
    }
}

const double Domain::WEIGHTSD2Q5   [ 5] = { 2./6., 1./6., 1./6., 1./6., 1./6 };
const double Domain::WEIGHTSD2Q9   [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Domain::WEIGHTSD3Q15  [15] = { 2./9., 1./9., 1./9., 1./9., 1./9.,  1./9.,  1./9., 1./72., 1./72. , 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
const double Domain::WEIGHTSD3Q19  [19] = { 1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
const size_t Domain::OPPOSITED2Q5  [ 5] = { 0, 3, 4, 1, 2 };                                                       ///< Opposite directions (D2Q5) 
const size_t Domain::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };                                           ///< Opposite directions (D2Q9) 
const size_t Domain::OPPOSITED3Q15 [15] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};                     ///< Opposite directions (D3Q15)
const size_t Domain::OPPOSITED3Q19 [19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};     ///< Opposite directions (D3Q19)
const Vec3_t Domain::LVELOCD2Q5  [ 5] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0} };
const Vec3_t Domain::LVELOCD2Q9  [ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Vec3_t Domain::LVELOCD3Q15 [15] =
{
	{ 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
	{ 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
	{-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} 
};

const Vec3_t Domain::LVELOCD3Q19 [19] =
{
	{ 0, 0, 0}, 
    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
    { 1, 1, 0}, {-1,-1, 0}, { 1, 0, 1}, {-1, 0,-1}, { 0, 1, 1}, { 0,-1,-1},
    { 1,-1, 0}, {-1, 1, 0}, { 1, 0,-1}, {-1, 0, 1}, { 0, 1,-1}, { 0,-1, 1}
};

const double   Domain::MD2Q9 [ 9][ 9] = { { 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0},
                                          {-4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0},{ 4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0},
                                          { 0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0},
                                          { 0.0, -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0},
                                          { 0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0},
                                          { 0.0,  0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0},
                                          { 0.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
                                          { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}};

const double Domain::MD3Q15 [15][15]  =  /*{ { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                            {-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                            {16.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                            { 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                            { 0.0,-4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                            { 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                            { 0.0, 0.0, 0.0,-4.0, 4.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 4.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                            { 0.0, 2.0, 2.0,-1.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0} };*/
                                        { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          {-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          {16.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          { 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                          { 0.0,-4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0,-4.0, 4.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 4.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                          { 0.0, 2.0, 2.0,-1.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                          { 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0} };

const double Domain::MD3Q19 [19][19]= { {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                        {-30.0, -11.0, -11.0, -11.0, -11.0, -11.0, -11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0},
                                        {12.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                        {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0},
                                        {0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0},
                                        {0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0},
                                        {0.0, 2.0, 2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0},
                                        {0.0, -4.0, -4.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0},
                                        {0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, -2.0, -2.0, 2.0, 2.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, -1.0, 1.0, 1.0, -1.0}};

}
#endif
