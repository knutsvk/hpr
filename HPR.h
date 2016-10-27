#ifndef __HPR_H
#define __HPR_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <Eigen/Core> // Matrix3d
#include <Eigen/Dense> // determinant() 
#include <boost/numeric/odeint.hpp> // ode integrator
#include <omp.h> // openmp parallelisation
#include <libconfig.h++> // configuration files

#include "SimpleArray.h"

/* TODO
 * vorticity 
 * stream function
 * initial conditions
 * - varying as f(x,y) in general
 * no-slip BC (reflective with const velocity)
 *  - laminar boundary layer
 *  - lid driven cavity
 * rethink naming of classes 
 * remove rubbish
 * comment code
 * add non-newtonian
 */

// Enumerators

enum Direction{ horizontal, vertical, radial };
enum BoundaryCondition{ transmissive, reflective, periodic, constant };

// Base class for Hyperbolic Peshkov Romenski based solver

class HyperbolicPeshkovRomenski
{
    protected: 
        double c_s; // shear wave speed
        double rho_0; // reference density 

        int nCellsX; // amount of cells in the x-direction
        int nCellsY; // amount of cells in the y-direction
        int nGhostCells; // ghost cells outside domain
        int nCellsTot; // total amount of cells
        double domain[4]; // xmin = domain[0], xmax = domain[1], ymin = domain[2],
        double dx; // cell width in x-direction
        double dy; // cell width in y-direction

        std::vector< SimpleArray< double, 14 > > consVars; // conserved variables 

        void xFlux( const SimpleArray< double, 14 >& Q, 
                SimpleArray< double, 14 >& F );
        void yFlux( const SimpleArray< double, 14 >& Q, 
                SimpleArray< double, 14 >& G );
        void forceFlux( double dt, double dr, int dir, 
                const SimpleArray< double, 14 >& Q_L, 
                const SimpleArray< double, 14 >& Q_R, 
                SimpleArray< double, 14 >& F );
        void slicFlux ( double dt, double dr, int dir, 
                const SimpleArray< double, 14 >& Q_2L, 
                const SimpleArray< double, 14 >& Q_L, 
                const SimpleArray< double, 14 >& Q_R, 
                const SimpleArray< double, 14 >& Q_2R, 
                SimpleArray< double, 14 >& F );
        void nonconservativeTerms( double dt, double dr, int dir, 
                const SimpleArray< double, 14 >& Q_L, 
                const SimpleArray< double, 14 >& Q_0, 
                const SimpleArray< double, 14 >& Q_R, 
                SimpleArray< double, 14 >& N );

    public:
        HyperbolicPeshkovRomenski( 
                double _shearSoundSpeed, double _referenceDensity, 
                int _nCellsX, int _nCellsY, double _domain[4] );
        virtual ~HyperbolicPeshkovRomenski(){}

        double getDensity( const SimpleArray< double, 14 >& Q );
        SimpleArray< double, 3 > getVelocity( const SimpleArray< double, 14 >& Q );
        Eigen::Matrix3d getDistortion( const SimpleArray< double, 14 >& Q );
        double getEnergy( const SimpleArray< double, 14 >& Q );

        virtual double getPressure( const SimpleArray< double, 14 >& Q ) = 0;
        Eigen::Matrix3d getShearStress( const SimpleArray< double, 14>& Q );

        virtual double microEnergy( double density, double pressure ) = 0;
        double mesoEnergy( Eigen::Matrix3d distortion );
        double macroEnergy( SimpleArray< double, 3 > velocity );

        virtual double getTimeStep( double c_CFL ) = 0;

        void initialise( double initDiscontPos, Direction initDiscontDir,
                double density[2], SimpleArray< double, 3 > velocity[2],
                Eigen::Matrix3d distortion[2], double pressure[2] ); 
                
        void boundaryConditions( BoundaryCondition type[4] );
        void xSweep( double dt );
        void ySweep( double dt );
        void diffuse();
        void renormalizeDistortion();
        void output2D( char* filename );
        void output1DSliceX( char* filename );
        void output1DSliceY( char* filename );

        bool isPhysical();

        // Functions taht will probably be deleted/replaced at some point: 
        Eigen::Matrix3d getCurlTerm( const SimpleArray< double, 14 >& Q_0, 
                const SimpleArray< double, 14 >& Q_L, 
                const SimpleArray< double, 14 >& Q_R,
                const SimpleArray< double, 14 >& Q_B, 
                const SimpleArray< double, 14 >& Q_T );
        void printDomain();
        void periodicBoundaryConditions();
        void initialiseDoubleShearLayer();
        void initialiseConvergenceTest();
        void exactConvergenceSolution();
};

// Fluid class, inherits HPR

class HPR_Fluid: public HyperbolicPeshkovRomenski
{
    protected:
        double gamma; // adiabatic constant = c_V / c_p
        double tau; // strain dissipation time / particle settled lifetime
        
    public:
        HPR_Fluid( double _shearSoundSpeed, double _referenceDensity, 
                int _nCellsX, int _nCellsY, double _domain[4], 
                double _gamma, double _strainDissipationTime );
        virtual ~HPR_Fluid(){}

        virtual double getPressure( const SimpleArray< double, 14 >& Q );
        Eigen::Matrix3d getSource( const SimpleArray< double, 14 >& Q );
        virtual double microEnergy( double density, double pressure );

        virtual double getTimeStep( double c_CFL );

        void integrateODE( double dt );
};

// Solid class, inherits HPR

class HPR_Solid: public HyperbolicPeshkovRomenski
{
    private: 
        double c_0;     // adiabatic sound speed
        double Gamma_0; // Grüneisen parameter
        double s_H;  // linear Hugoniot slope coefficient

    public:
        HPR_Solid( double _shearSoundSpeed, double _referenceDensity, 
                int _nCellsX, int _nCellsY, double _domain[4], 
                double _c_0, double _Gamma_0, double _s_H );
        virtual ~HPR_Solid(){}

        virtual double getPressure( const SimpleArray< double, 14 >& Q );
        virtual double microEnergy( double density, double pressure );

        virtual double getTimeStep( double c_CFL );
};

// ODE struct
typedef SimpleArray< double, 14 > state_type;
typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > stepper_type;

struct HPR_ODE
{
    double tau;
    double rho_0;
    HPR_ODE(double _tau, double rho_0);
    void operator()( const state_type& Q, state_type& S, const double t );
};

// Functions independent of classes

double slopeLimiter( double q_minus, double q_0, double q_plus );
void configurate( const char* inputFile, int& nCellsX, int& nCellsY,
        double& CFL, double& tStop, double& c_s, double& rho_0, double& gamma,
        double& tau, double domain[4], double& initDiscontPos,
        Direction& initDiscontDir, BoundaryCondition BCs[4], double
        initDensity[2], SimpleArray< double, 3 > initVelocity[2],
        Eigen::Matrix3d initDistortion[2], double initPressure[2] );

#endif
