#ifndef __HPR_H
#define __HPR_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <omp.h> // openmp parallelisation
#include <libconfig.h++> // configuration files

#include "SimpleArray.h"

// Enumerators

enum Direction{ horizontal, vertical, radial };
enum BoundaryCondition{ transmissive, reflective, periodic, constant };

// Base class for 2D solver of Euler equations

class Euler
{
    private: 
        int nCellsX;        // amount of cells in the x-direction
        int nCellsY;        // amount of cells in the y-direction
        int nGhostCells;    // ghost cells outside domain
        int nCellsTot;      // total amount of cells

        double domain[4];   // xmin = domain[0], xmax = domain[1], ymin = domain[2],
        double dx;          // cell width in x-direction
        double dy;          // cell width in y-direction
        double gamma;       // adiabatic constant = c_V / c_p

        std::vector< SimpleArray< double, 4 > > consVars; // conserved variables 

        void xFlux( const SimpleArray< double, 4 >& Q, 
                SimpleArray< double, 4 >& F );
        void yFlux( const SimpleArray< double, 4 >& Q, 
                SimpleArray< double, 4 >& G );
        void forceFlux( double dt, double dr, int dir, 
                const SimpleArray< double, 4 >& Q_L, 
                const SimpleArray< double, 4 >& Q_R, 
                SimpleArray< double, 4 >& F );
        void slicFlux ( double dt, double dr, int dir, 
                const SimpleArray< double, 4 >& Q_2L, 
                const SimpleArray< double, 4 >& Q_L, 
                const SimpleArray< double, 4 >& Q_R, 
                const SimpleArray< double, 4 >& Q_2R, 
                SimpleArray< double, 4 >& F );

    public:
        Euler( int _nCellsX, int _nCellsY, double _domain[4], double _gamma );
        virtual ~Euler(){}

        double getDensity( const SimpleArray< double, 4 >& Q );
        SimpleArray< double, 2 > getVelocity( const SimpleArray< double, 4 >& Q );
        double getEnergy( const SimpleArray< double, 4 >& Q );
        double getPressure( const SimpleArray< double, 4 >& Q );
        double microEnergy( double density, double pressure );
        double macroEnergy( SimpleArray< double, 2 > velocity );

        double getTimeStep( double c_CFL );

        void initialise( double initDiscontPos, Direction initDiscontDir,
                double density[2], SimpleArray< double, 2 > velocity[2],
                double pressure[2] ); 
                
        void boundaryConditions( BoundaryCondition type[4] );
        void xSweep( double dt );
        void ySweep( double dt );

        void output2D( char* filename );
        void output1DSliceX( char* filename );
        void output1DSliceY( char* filename );

        bool isPhysical();

        // Functions taht will probably be deleted/replaced at some point: 
        void initialiseConvergenceTest();
        void exactConvergenceSolution();
};

// Functions independent of classes

double slopeLimiter( double q_minus, double q_0, double q_plus );
void configurate( const char* inputFile, int& nCellsX, int& nCellsY, 
        double& CFL, double& tStop, double& gamma, double domain[4], 
        double& initDiscontPos, Direction& initDiscontDir, 
        BoundaryCondition BCs[4], double initDensity[2], 
        SimpleArray< double, 2 > initVelocity[2], double initPressure[2] );

#endif
