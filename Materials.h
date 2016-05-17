#ifndef __MATERIALS_H
#define __MATERIALS_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

#include "SimpleArray.h"

/* TODO
 * make classes which inherit Material: Fluid and Solid
 * allow for sources
 * extend to 2D
 * add non-conservative
 * try viscous
 */

class Material
{
    private: 
        unsigned nCells;
        const static unsigned nGhostCells = 2;
        double domain[2];
        double dx;
        double gamma;

    protected: 
        std::vector< SimpleArray< double, 3 > > consVars;
        std::vector< SimpleArray< double, 3 > > xDirFlux;

        void flux(const SimpleArray< double, 3 >& Q, SimpleArray< double, 3 >& F);

    public:
        Material(const int _nCells, const double _domain[2]);
        virtual ~Material(){}

        std::vector<double> getDensity();
        std::vector<double> getVelocity();
        std::vector<double> getPressure();
        std::vector<double> getInternalEnergy();
        double timeStep(const double c_CFL);

        void initialize(const double interfacePos, const double density[2],
                const double velocity[2], const double pressure[2]);
        void transmissiveBCs();
        void reflectiveBCs();
        void force(double dt);
        void slic(double dt);
        void advancePDE(const double dt);
        void output();
};

double slopeLimiter(double q_minus, double q_0, double q_plus);

#endif
