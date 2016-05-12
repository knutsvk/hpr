#ifndef __MATERIALS_H
#define __MATERIALS_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <valarray>

class Material
{
    private: 
        int nCells;
        int nGhostCells;
        double domain[2];
        double dx;
        double gamma;       // TODO:Make specific to idealGas

    protected: 
        // TODO: test kvector and blitz++
        std::vector<std::valarray<double> > primVars;
        std::vector<std::valarray<double> > consVars;
        std::vector<std::valarray<double> > tempVars;
        std::vector<std::valarray<double> > xDirFlux;

        void updatePrimitive();
        void updateConserved();

        void flux(const std::valarray<double>& Q, std::valarray<double>& F);

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
