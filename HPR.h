#ifndef __HPR_H
#define __HPR_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <Eigen/Core> // Matrix3d
#include <Eigen/Dense> // determinant() 
#include <boost/numeric/odeint.hpp> // ode integrator
#include <omp.h> // openmp parallelisation

#include "SimpleArray.h"

/* TODO
 * omp for subcycling in IntegrateODE()
 * change state_type in odeint
 * add non-conservative
 * extend to 2D, test with cylindrical explosion
 * viscous test cases
 *  - first problem of stokes
 *  - laminar boundary layer
 *  - lid driven cavity
 */

class HyperbolicPeshkovRomenski
{
    protected: 
        double c_s; // shear wave speed
        double rho_0; // reference density 

        unsigned nCells; // amount of cells in the computational domain
        const static unsigned nGhostCells = 2; // ghost cells outside computational domain
        double domain[2]; // xmin = domain[0], xmax = domain[1]
        double dx; // cell width

        std::vector< SimpleArray< double, 14 > > consVars; // conserved variables in each cell
        std::vector< SimpleArray< double, 14 > > xDirFlux; // fluxes at each cell boundary 

        void flux( const SimpleArray< double, 14 >& Q, 
                SimpleArray< double, 14 >& F );

    public:
        HyperbolicPeshkovRomenski( 
                double _shearSoundSpeed, double _referenceDensity, 
                int _nCells, double _domain[2] );
        virtual ~HyperbolicPeshkovRomenski(){}

        double getDensity( const SimpleArray< double, 14 >& Q );
        // TODO: use std::vector or std::array (or double[3]) instead? 
        SimpleArray< double, 3 > getVelocity( const SimpleArray< double, 14 >& Q );
        Eigen::Matrix3d getDistortion( const SimpleArray< double, 14 >& Q );
        double getEnergy( const SimpleArray< double, 14 >& Q );

        virtual double getPressure( const SimpleArray< double, 14 >& Q ) = 0;
        Eigen::Matrix3d getShearStress( const SimpleArray< double, 14>& Q );

        virtual double microEnergy( double density, double pressure ) = 0;
        double mesoEnergy( Eigen::Matrix3d distortion );
        double macroEnergy( SimpleArray< double, 3 > velocity );

        virtual double getTimeStep( double c_CFL ) = 0;

        void initialize( double initDiscontPos, double density_L, double density_R, 
                SimpleArray< double, 3 > velocity_L, SimpleArray< double, 3 > velocity_R, 
                Eigen::Matrix3d distortion_L, Eigen::Matrix3d distortion_R, 
                double pressure_L, double pressure_R ); 
                
        void transmissiveBCs();
        void reflectiveBCs();
        void force( double dt );
        void slic( double dt );
        void advancePDE( double dt );
        void renormalizeDistortion();
        void output();
};

class HPR_Fluid: public HyperbolicPeshkovRomenski
{
    protected:
        double gamma; // adiabatic constant = c_V / c_p
        double tau; // strain dissipation time, alternatively particle settled lifetime
        
    public:
        HPR_Fluid( double _shearSoundSpeed, double _referenceDensity, 
                int _nCells, double _domain[2], 
                double _gamma, double _strainDissipationTime );
        virtual ~HPR_Fluid(){}

        virtual double getPressure( const SimpleArray< double, 14 >& Q );
        Eigen::Matrix3d getSource( const SimpleArray< double, 14 >& Q );
        virtual double microEnergy( double density, double pressure );

        virtual double getTimeStep( double c_CFL );

        void integrateODE( double dt );
};

class HPR_Solid: public HyperbolicPeshkovRomenski
{
    private: 
        double c_0;     // adiabatic sound speed
        double Gamma_0; // Grüneisen parameter
        double s_H;  // linear Hugoniot slope coefficient

    public:
        HPR_Solid( double _shearSoundSpeed, double _referenceDensity, 
                int _nCells, double _domain[2], 
                double _c_0, double _Gamma_0, double _s_H );
        virtual ~HPR_Solid(){}

        virtual double getPressure( const SimpleArray< double, 14 >& Q );
        virtual double microEnergy( double density, double pressure );

        virtual double getTimeStep( double c_CFL );
};

// ODE struct
// TODO: Use SimpleArray< double, 14 > as state_type
typedef std::vector< double > state_type; 
typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > stepper_type;

struct HPR_ODE
{
    double tau;
    double rho_0;
    HPR_ODE(double _tau, double rho_0);
    void operator()( const state_type& Q, state_type& S, const double t );
};

double slopeLimiter( double q_minus, double q_0, double q_plus );

#endif
