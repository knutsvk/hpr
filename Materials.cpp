#ifndef __MATERIALS_CPP
#define __MATERIALS_CPP

#include "Materials.h"

void Material::flux(const SimpleArray< double, 3 >& Q, SimpleArray< double, 3 >& F)
{
    F[0] = Q[1];
    F[1] = (gamma - 1.0) * Q[2] + 0.5 * (3.0 - gamma) * Q[1] * Q[1] / Q[0];
    F[2] = (gamma * Q[2] - 0.5 * (gamma - 1.0) * Q[1] * Q[1] / Q[0]) * Q[1] / Q[0];
}

Material::Material(const int _nCells, const double _domain[2])
{
    for(unsigned i = 0; i < 2; i++)
    {
        domain[i] = _domain[i];
    }

    dx = (domain[1] - domain[0]) / nCells;

    consVars.resize(_nCells + 2 * nGhostCells);
    xDirFlux.resize(_nCells + 1);
}

std::vector<double> Material::getDensity()
{
    std::vector<double> rho(nCells + 2 * nGhostCells);
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        rho[i] = consVars[i][0];
    }
    return rho;
}

std::vector<double> Material::getVelocity()
{
    std::vector<double> u(nCells + 2 * nGhostCells);
    for(unsigned i = nGhostCells; i < nCells + 2 * nGhostCells; i++)
    {
        u[i] = consVars[i][1] / consVars[i][0];
    }
    return u;
}

std::vector<double> Material::getPressure()
{
    std::vector<double> p(nCells + 2 * nGhostCells);
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        p[i] = (consVars[i][2] - 0.5 * consVars[i][1] *
                consVars[i][1] / consVars[i][0]) * (gamma - 1.0);
    }
    return p;
}

std::vector<double> Material::getInternalEnergy()
{
    std::vector<double> e(nCells + 2 * nGhostCells);
    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        e[i] = (consVars[i][2] - 0.5 * consVars[i][1] *
                consVars[i][1] / consVars[i][0]) / (consVars[i][0]);
    }
    return e;
}

double Material::timeStep(const double c_CFL)
{
    std::vector<double> rho = getDensity();
    std::vector<double> u = getVelocity();
    std::vector<double> p = getPressure();
    std::vector<double> S(nCells + 2 * nGhostCells); 

    for(unsigned i = 0; i < nCells + 2 * nGhostCells; i++)
    {
        S[i] = fabs(u[i]) 
            + sqrt(gamma * p[i] / rho[i]);
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
            consVars[i][0] = density[0];
            consVars[i][1] = density[0] * velocity[0];
            consVars[i][2] = 0.5 * density[0] * velocity[0] * velocity[0] +
                pressure[0] / (gamma - 1.0);
        }
        else
        {
            consVars[i][0] = density[1];
            consVars[i][1] = density[1] * velocity[1];
            consVars[i][2] = 0.5 * density[1] * velocity[1] * velocity[1] +
                pressure[1] / (gamma - 1.0);
        }
    }
}

void Material::transmissiveBCs()
{
    for(unsigned i = 0; i < nGhostCells; i++)
    {
        consVars[i] = consVars[i + 1];
        consVars[nCells + 2 * nGhostCells -1 - i] = 
            consVars[nCells + 2 * nGhostCells - 1 - (i + 1)];
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

void Material::force(double dt)
{
    int L, R;
    SimpleArray< double, 3> F_L, F_R, Q_0, F_0;
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

void Material::advancePDE(const double dt)
{
    std::vector< SimpleArray< double, 3 > > tempVars(nCells + 2 * nGhostCells);
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
    std::vector<double> rho, u, p, e;
    rho = getDensity();
    u = getVelocity();
    p = getPressure();
    e = getInternalEnergy();

    std::cout << "x \trho \tu \tp \te" << std::endl;
    for(unsigned i = nGhostCells; i < nCells +  nGhostCells; i++)
    {
        std::cout << domain[0] + (i - nGhostCells + 0.5) * dx << " \t" 
            << rho[i] << " \t" << u[i] << " \t" << p[i] << " \t" << e[i] 
            << std::endl;
    }
}

#endif
