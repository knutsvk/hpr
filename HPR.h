#ifndef __HPR_H
#define __HPR_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <omp.h>

#include "SimpleArray.h"

/* TODO
 * make classes which inherit Material: Fluid and Solid
 * allow for sources
 * extend to 2D
 * add non-conservative
 * try viscous
 */

class HyperbolicPeshkovRomenski
{
    protected: 
        double c_s; 
        double tau; 
        double rho_0;

        unsigned nCells;
        const static unsigned nGhostCells = 2;
        double domain[2];
        double dx;

        std::vector< SimpleArray< double, 14 > > consVars;
        std::vector< SimpleArray< double, 14 > > xDirFlux;

        void flux( const SimpleArray< double, 14 >& Q, 
                SimpleArray< double, 14 >& F );

    public:
        HyperbolicPeshkovRomenski( double _shearSoundSpeed, 
                double _strainDissipationTime, double _referenceDensity, 
                int _nCells, double _domain[2] );
        virtual ~HyperbolicPeshkovRomenski(){}

        double getDensity( const SimpleArray< double, 14 >& Q );
        double* getVelocity( const SimpleArray< double, 14 >& Q );
        Eigen::Matrix3d getDistortion( const SimpleArray< double, 14 >& Q );
        double getEnergy( const SimpleArray< double, 14 >& Q );

        virtual double getPressure( const SimpleArray< double, 14 >& Q ) = 0;
        Eigen::Matrix3d getShearStress( const SimpleArray< double, 14>& Q );
        Eigen::Matrix3d getSource( const SimpleArray< double, 14 >& Q );

        virtual double microEnergy( double density, double pressure ) = 0;
        double mesoEnergy( Matrix3d distortion );
        double macroEnergy( double* velocity );

        double getTimeStep( double c_CFL );

        void initialize( double initDiscontPos, double density_L, double density_R, 
                Eigen::Vector3d velocity_L, Eigen::Vector3d velocity_R, 
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
        double gamma;
        
    public:
        HPR_Fluid( double _shearSoundSpeed, double _strainDissipationTime,
                double _referenceDensity, int _nCells, double _domain[2],
                double _gamma );
        virtual ~HPR_Fluid(){}

        virtual double getPressure( const SimpleArray< double, 14 >& Q );
        virtual double microEnergy( double density, double pressure );

};

class HPR_Solid: public HyperbolicPeshkovRomenski
{
    private: 
        double c_0; 
        double Gamma_0;

    public:
        HPR_Solid( double _shearSoundSpeed, double _strainDissipationTime,
                double _referenceDensity, int _nCells, double _domain[2],
                double _rho_0, double _c_0, double _Gamma_0 );
        virtual ~HPR_Solid(){}

        virtual double getPressure( const SimpleArray< double, 14 >& Q );
        virtual double microEnergy( double density, double pressure );
};

double slopeLimiter( double q_minus, double q_0, double q_plus );

#endif
