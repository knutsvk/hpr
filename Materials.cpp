#ifndef __MATERIALS_CPP
#define __MATERIALS_CPP

#include "Materials.h"

void Material::updatePrimitive()
{
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        primVars[i][0] = consVars[i][0];
        primVars[i][1] = consVars[i][1] / consVars[i][0];
        primVars[i][2] = (gamma - 1.0) * (consVars[i][2] - 0.5 * consVars[i][1]
                * consVars[i][1] / consVars[i][0]);
    }
}

void Material::updateConserved()
{
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        consVars[i][0] = primVars[i][0];
        consVars[i][1] = primVars[i][0] * primVars[i][1];
        consVars[i][2] = 0.5 * primVars[i][0] * primVars[i][1] * primVars[i][1]
            + primVars[i][2] / (gamma - 1.0);
    }
}

void Material::flux(const std::valarray<double>& Q, std::valarray<double>& F)
{
    F[0] = Q[1];
    F[1] = (gamma - 1.0) * Q[2] + 0.5 * (3.0 - gamma) * Q[1] * Q[1] / Q[0];
    F[2] = (gamma * Q[2] - 0.5 * (gamma - 1.0) * Q[1] * Q[1] / Q[0]) * Q[1] / Q[0];
}

void Material::force(double dt)
{
    int L, R;
    std::valarray<double> F_L(3), F_R(3), Q_0(3), F_0(3);
    for(int i = 0; i < nCells + 1; i++)
    {
        L = i + nGhostCells - 1;
        R = i + nGhostCells;
        flux(consVars[L], F_L);
        flux(consVars[R], F_R);
        Q_0 = 0.5 * (consVars[L] + consVars[R]) 
            + 0.5 * dt / dx * (F_L - F_R);
        flux(Q_0, F_0);
        xDirFlux[i] = 0.5 * (F_0 + 0.5 * (F_L + F_R)) 
            + 0.25 * dx / dt * (consVars[L] - consVars[R]); 
    }
}

Material::Material(const int _nCells, const double _domain[2])
{
    nCells = _nCells;
    nGhostCells = 2;

    for(unsigned i = 0; i < 2; i++)
    {
        domain[i] = _domain[i];
    }

    dx = (domain[1] - domain[0]) / nCells;

    primVars.resize(nCells + 2 * nGhostCells);
    consVars.resize(nCells + 2 * nGhostCells);
    tempVars.resize(nCells + 2 * nGhostCells);
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        primVars[i].resize(3);
        consVars[i].resize(3);
        tempVars[i].resize(3);
    }

    xDirFlux.resize(nCells + 1);
    for(unsigned i = 0; i < nCells + 1; i++)
    {
        xDirFlux[i].resize(3);
    }

    gamma = 1.4;
}

double Material::timeStep(const double c_CFL)
{
    std::vector<double> S(nCells + 2 * nGhostCells); 
    updatePrimitive();

    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        S[i] = std::fabs(primVars[i][1]) 
            + sqrt(gamma * primVars[i][2] / primVars[i][0]);
    }

    return c_CFL * dx / *std::max_element(S.begin(), S.end());
}

void Material::initialize(const double initDiscontPos, const double density[2],
        const double velocity[2], const double pressure[2])
{
    double x;

    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        x = (i - nGhostCells + 0.5) * dx;

        if(x < initDiscontPos)
        {
            primVars[i][0] = density[0];
            primVars[i][1] = velocity[0];
            primVars[i][2] = pressure[0];
        }
        else
        {
            primVars[i][0] = density[1];
            primVars[i][1] = velocity[1];
            primVars[i][2] = pressure[1];
        }
    }
    updateConserved();
}

void Material::transmissiveBCs()
{
    for(unsigned i = 0; i < nGhostCells; i++)
    {
        consVars[i] = consVars[i + 1];
        consVars[nCells + 2 * nGhostCells - i] = 
            consVars[nCells + 2 * nGhostCells - (i + 1)];
    }
}

void Material::reflectiveBCs()
{
    for(unsigned i = 0; i < nGhostCells; i++)
    {
        consVars[i] = consVars[i + 1];
        consVars[i][1] *= -1.0;
        consVars[nCells + 2 * nGhostCells - i] = 
            consVars[nCells + 2 * nGhostCells - (i + 1)];
        consVars[nCells + 2 * nGhostCells - i][1] *= -1;
    }
}

void Material::advancePDE(const double dt)
{
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        tempVars[i] = consVars[i];
    } 
    
    int L, R;
    for(unsigned i = nGhostCells; i < nCells + nGhostCells; i++)
    {
        L = i - nGhostCells;
        R = i - nGhostCells + 1;
        consVars[i] = tempVars[i] + dt / dx * 
            (xDirFlux[L] - xDirFlux[R]);
    } 
}

void Material::output()
{
    updatePrimitive();

    std::cout << "x \trho \tu \tp \te" << std::endl;
    for(unsigned i = nGhostCells; i < nCells +  nGhostCells; i++)
    {
        std::cout << (i - nGhostCells + 0.5) * dx << " \t" << primVars[i][0] << " \t"
            << primVars[i][1] << " \t" << primVars[i][2] << " \t"
            << primVars[i][2] / ((gamma - 1.0) * primVars[i][0]) << std::endl;
    }
}

#endif
