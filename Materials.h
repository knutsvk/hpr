#ifndef __MATERIALS_H
#define __MATERIALS_H

#include <algorithm>
#include <iostream>
#include <vector>
#include "SimpleArray.h"

class Material
{
    private: 
        int nCells;
        int nGhostCells;
        double domain[2];
        double dx;
        double gamma;       // TODO:Make specific to idealGas

    protected: 
        std::vector< SimpleArray< double, 3 > > primVars;
        std::vector< SimpleArray< double, 3 > > consVars;
        std::vector< SimpleArray< double, 3 > > tempVars;
        std::vector< SimpleArray< double, 3 > > xDirFlux;

        void updatePrimitive();
        void updateConserved();

        void flux(const SimpleArray< double, 3 >& Q, SimpleArray< double, 3 >& F);

    public:
        Material(const int _nCells, const double _domain[2]);
        virtual ~Material(){}

        double timeStep(const double c_CFL);

        void initialize(const double interfacePos, const double density[2],
                const double velocity[2], const double pressure[2]);
        void transmissiveBCs();
        void reflectiveBCs();
        void force(double dt);
        void advancePDE(const double dt);
        void output();
};

#endif
