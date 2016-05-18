#ifndef __HPR_CPP
#define __HPR_CPP

#include "HPR.h"

void HyperbolicPeshkovRomenski::flux( 
        const SimpleArray< double, 14 >& Q, 
        SimpleArray< double, 14 >& F )
{
    double rho = getDensity( Q );
    double u [3] = getVelocity( Q );
    Matrix3d A = getDistortion( Q );
    double E = getEnergy( Q );

    Matrix3d sigma = shearStress( Q );
    double p = pressure( Q );
    
    F[0] = rho * u[0];
    F[1] = rho * u[0] * u[0] - sigma(0, 0) + p;
    F[2] = rho * u[0] * u[1] - sigma(1, 0);
    F[3] = rho * u[0] * u[2] - sigma(2, 0);
    F[4] = A( 0, 0 ) * u[0] + A(0, 1) * u[1] + A(0, 2) * u [2];
    F[5] = 0.0;
    F[6] = 0.0;
    F[7] = A(1, 0) * u[0] + A(1, 1) * u[1] + A(1, 2) * u [2];
    F[8] = 0.0;
    F[9] = 0.0;
    F[10] = A(2, 0) * u[0] + A(2, 1) * u[1] + A(2, 2) * u [2];
    F[11] = 0.0;
    F[12] = 0.0;
    F[13] = rho * u[0] * E + u[0] * p - sigma(0, 0) * u[0] 
        - sigma(1, 0) * u[1] - sigma(2, 0) * u[2];
}

HyperbolicPeshkovRomenski::HyperbolicPeshkovRomenski( 
        double _shearSoundSpeed, double _strainDissipationTime, 
        double _referenceDensity, int _nCells, double _domain[2] )
{
    c_s = _shearSoundSpeed; 
    tau = _strainDissipationTime; 
    rho_0 = _referenceDensity;
    nCells = _nCells;

    for(unsigned i = 0; i < 2; i++)
    {
        domain[i] = _domain[i];
    }

    dx = ( domain[1] - domain[0] ) / nCells;

    consVars.resize( nCells + 2 * nGhostCells );
    xDirFlux.resize( nCells + 1 );
}

double HyperbolicPeshkovRomenski::getDensity( 
        const SimpleArray< double, 14 >& Q )
{
    return Q[0];
}

double* HyperbolicPeshkovRomenski::getVelocity( 
        const SimpleArray< double, 14 >& Q )
{
    double u[3]; 
    for( unsigned i = 0; i < 3; i++ ) 
    { 
        u[i] = Q[i + 1] / Q[0];
    }
    return u[0];
}

Eigen::Matrix3d HyperbolicPeshkovRomenski::getDistortion( 
        const SimpleArray< double, 14 >& Q )
{
    Matrix3d A;
    for( unsigned i = 0; i < 3; i++ )
    {
        for( unsigned j = 0; j < 3; j++ )
        {
            A(i, j) = Q[4 + 3 * i + j];
        }
    }
    return A;
}

double HyperbolicPeshkovRomenski::getEnergy(
        const SimpleArray< double, 14 >& Q )
{
    return Q[13] / Q[0];
}

Eigen::Matrix3d HyperbolicPeshkovRomenski::getShearStress( 
        const SimpleArray< double, 14 >& Q )
{
    double rho = getDensity( Q );
    Eigen::Matrix3d A = getDistortion( Q );
    Eigen::Matrix3d G = A.transpose() * A;
    Eigen::Matrix3d devG = G - G.trace() * Eigen::Matrix3d::Identity() / 3.0;
    return - rho * c_s * c_s * G * devG;
}

Eigen::Matrix3d HyperbolicPeshkovRomenski::getSource(
        const SimpleArray< double, 14 >& Q )
{
    double nu = getDensity( Q ) / rho_0;
    Eigen::Matrix3d A = getDistortion( Q );
    Eigen::Matrix3d G = A.transpose() * A;
    Eigen::Matrix3d devG = G - G.trace() * Eigen::Matrix3d::Identity() / 3.0;
    return - 3.0 * pow(nu, 5.0 / 3.0) * A * devG / tau;
}

double mesoEnergy( Matrix3d distortion )
{
    Matrix3d G = distortion.trace() * distortion;
    return 0.25 * c_s * c_s * ( G.trace() - ( G * G ).trace() / 3.0 );
}

double macroEnergy( double* velocity )
{
    return 0.5 * ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
}

double HyperbolicPeshkovRomenski::timeStep( const double c_CFL )
{
    std::vector<double> S( nCells + 2 * nGhostCells ); 
    double rho;
    double* u; 
    double p;

#pragma omp parallel for private( rho, u, p )
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        rho = getDensity( consVars[i] );
        u = getVelocity( consVars[i] );
        p = getPressure( consVars[i] );
        S[i] = fabs( u[0] ) + sqrt( gamma * p / rho + 4.0 * c_s * c_s / 3.0 );
    }

    return c_CFL * dx / *std::max_element( S.begin(), S.end() );
}

void HyperbolicPeshkovRomenski::initialize( double initDiscontPos, 
        double density_L, double density_R, 
        double* velocity_L, double* velocity_R, 
        Eigen::Matrix3d distortion_L, Eigen::Matrix3d distortion_R,
        double pressure_L, double pressure_R )
{
    double x;
#pragma omp parallel for private( x )
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        x = domain[0] + ( i - nGhostCells + 0.5 ) * dx;

        if(x < initDiscontPos)
        {
            consVars[i][0] = density_L;
            for( unsigned j = 0; j < 3; j++ )
            {
                consVars[i][j + 1] = density_L * velocity_L[j];
                for( unsigned k = 0; k < 3; k++ )
                {
                    consVars[i][4 + 3 * j + k] = distortion_L(j, k);
                }
            }
            consVars[i][13] = rho_L * 
                ( microEnergy( density_L, pressure_L ) +
                  mesoEnergy( distortion_L ) + 
                  macroEnergy( velocity_L ) );
        }
        else
        {
            consVars[i][0] = density_R;
            for( unsigned j = 0; j < 3; j++ )
            {
                consVars[i][j + 1] = density_R * velocity_R[j];
                for( unsigned k = 0; k < 3; k++ )
                {
                    consVars[i][4 + 3 * j + k] = distortion_R(j, k);
                }
            }
            consVars[i][13] = rho_R * 
                ( microEnergy( density_R, pressure_R ) +
                  mesoEnergy( distortion_R ) + 
                  macroEnergy( velocity_R ) );
        }
    }
}

void HyperbolicPeshkovRomenski::transmissiveBCs()
{
    for( unsigned i = 0; i < nGhostCells; i++ )
    {
        consVars[i] = consVars[i + 1];
        consVars[nCells + 2 * nGhostCells - 1 - i] = 
            consVars[nCells + 2 * nGhostCells - 1 - (i + 1)];
    }
}

void HyperbolicPeshkovRomenski::reflectiveBCs()
{
    for( unsigned i = 0; i < nGhostCells; i++ )
    {
        consVars[i] = consVars[i + 1];
        consVars[i][1] *= - 1.0;
        consVars[nCells + 2 * nGhostCells - i] = 
            consVars[nCells + 2 * nGhostCells - (i + 1)];
        consVars[nCells + 2 * nGhostCells - i][1] *= - 1.0;
    }
}

void HyperbolicPeshkovRomenski::force( double dt )
{
    int L, R;
    SimpleArray< double, 14 > F_L, F_R, Q_0, F_0;

#pragma omp parallel for private( L, R, F_L, F_R, Q_0, F_0 )
    for( unsigned i = 0; i < nCells + 1; i++ )
    {
        L = i + nGhostCells - 1;
        R = i + nGhostCells;

        flux( consVars[L], F_L );
        flux( consVars[R], F_R );

        Q_0 = 0.5 * ( consVars[L] + consVars[R] )
            + 0.5 * dt / dx * ( F_L - F_R );
        flux( Q_0, F_0 );

        xDirFlux[i] = 0.5 * ( F_0 + 0.5 * ( F_L + F_R ) ) 
            + 0.25 * dx / dt * ( consVars[L] - consVars[R] ); 
    }
}

void HyperbolicPeshkovRomenski::slic( double dt )
{
    int L, R;
    SimpleArray< double, 14 > xi_L, xi_R, 
        Q_L_plus, Q_R_plus, Q_L_0, Q_R_0,
        F_L_plus, F_R_plus, F_L_0, F_R_0, 
        Q_L, Q_R, Q_0, 
        F_L, F_R, F_0;

#pragma omp parallel for private( L, R, xi_L, xi_R, Q_L_plus, Q_R_plus, Q_L_0,\
        Q_R_0, F_L_plus, F_R_plus, F_L_0, F_R_0, Q_L, Q_R, Q_0, F_L, F_R, F_0 )

    for( unsigned i = 0; i < nCells + 1; i++ )
    {
        L = i + nGhostCells - 1;
        R = i + nGhostCells;

        // Calculate TVD slope limiters
        for( unsigned j = 0; j < 14; j++ )
        {
            xi_L[j] = slopeLimiter( consVars[L - 1][j], consVars[L][j],
                    consVars[R][j] );
            xi_R[j] = slopeLimiter( consVars[L][j], consVars[R][j], 
                    consVars[R + 1][j] );
        }

        // Boundary extrapolated values
        Q_L_plus = consVars[R] - 0.25 * xi_R * ( consVars[R + 1] - consVars[L] );
        Q_R_plus = consVars[R] + 0.25 * xi_R * ( consVars[R + 1] - consVars[L] );
        Q_L_0 = consVars[L] - 0.25 * xi_L * ( consVars[R] - consVars[L - 1] );
        Q_R_0 = consVars[L] + 0.25 * xi_L * ( consVars[R] - consVars[L - 1] );

        flux( Q_L_plus, F_L_plus );
        flux( Q_R_plus, F_R_plus );
        flux( Q_L_0, F_L_0 );
        flux( Q_R_0, F_R_0 );

        // Evolve by time 0.5 * dt
        Q_R = Q_L_plus + 0.5 * dt / dx * ( F_L_plus - F_R_plus );
        Q_L = Q_R_0 + 0.5 * dt / dx * ( F_L_0 - F_R_0 );

        // FORCE flux
        flux( Q_L, F_L );
        flux( Q_R, F_R );
        Q_0 = 0.5 * ( Q_L + Q_R ) 
            + 0.5 * dt / dx * ( F_L - F_R );
        flux( Q_0, F_0 );
        xDirFlux[i] = 0.5 * ( F_0 + 0.5 * ( F_L + F_R ) ) 
            + 0.25 * dx / dt * ( Q_L - Q_R ); 
    }
}

void HyperbolicPeshkovRomenski::advancePDE( double dt )
{
    std::vector< SimpleArray< double, 14 > > tempVars( nCells + 2 * nGhostCells );
#pragma omp parallel for 
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        tempVars[i] = consVars[i];
    } 
    
    int L, R;
    for( unsigned i = nGhostCells; i < nCells + nGhostCells; i++ )
    {
        L = i - nGhostCells;
        R = i - nGhostCells + 1;
        consVars[i] = tempVars[i] + dt / dx * 
            ( xDirFlux[L] - xDirFlux[R] );
    } 
}

void HyperbolicPeshkovRomenski::renormalizeDistortion()
{
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        double rho = getDensity( consVars[i] );
        Matrix3d A = getDistortion( consVars[i] );
        double scaleFactor = pow( rho / ( rho_0 * A.determinant() ), 1.0 / 3.0 );
        A *= scaleFactor;
        for( unsigned j = 0; j < 3; j++ )
        {
            for( unsigned k = 0; k < 3; k++ )
            {
                consVars[i][4 + 3 * j + k] = A(j, k);
            }
        }
    }
}

void HyperbolicPeshkovRomenski::output()
{
    double x; 
    double rho; 
    double u[3];
    double p; 
    Matrix3d sigma;

    // TODO: add entropy
    //
    std::cout << "x" << "\t" << "rho" << "\t" << "p" << "\t" 
        << "u" << "\t" << "v" << "\t" << "w" << "\t" 
        << "sigma11" << "\t" << "sigma12" << "\t" << "sigma13" << "\t" 
        << "sigma22" << "\t" << "sigma23" << "\t" << "sigma33" << std::endl; 

    for( unsigned i = nGhostCells; i < nCells +  nGhostCells; i++ )
    {
        rho = getDensity( consVars[i] );
        u = getVelocity( consVars[i] );
        p = getPressure( consVars[i] );
        sigma = getShearStress( consVars[i] );

        std::cout << domain[0] + ( i - nGhostCells + 0.5 ) * dx << " \t" 
            << rho << " \t" << p << " \t" << u[0] << " \t" << u[1] << " \t" 
            << u[2] << " \t" << sigma(0, 0) - p << "\t" << sigma(0, 1) << "\t"
            << sigma(0, 2) << "\t" << sigma(1, 1) - p << "\t" 
            << sigma(1, 2) << "\t" << sigma(2, 2) - p << "\t" << std::endl;
    }
}

HPR_Fluid 

double slopeLimiter(double q_min, double q_0, double q_plus)
{ // TODO: move to different file
    double Delta_L, Delta_R; 

    if(fabs(q_0 - q_min) < 1.0e-16) 
    {
        Delta_L = copysign(1.0e-16, q_0 - q_min);
    }
    else 
    {
        Delta_L = q_0 - q_min;
    }

    if(fabs(q_plus -q_0) < 1.0e-16) 
    {
        Delta_R = copysign(1.0e-16, q_plus - q_0);
    }
    else 
    {
        Delta_R = q_plus -q_0;
    }
    
    double r = Delta_L / Delta_R;

    // superbee
    if(r <= 0.0)
    {
        return 0.0;
    }
    else if(r <= 0.5) 
    {
        return 2.0 * r; 
    }
    else if(r <= 1.0) 
    {
        return 1.0; 
    }
    else 
    { //TODO: std::min_element()
        return (2.0 < r) ?  (2.0 < 2.0 / (1.0 + r) ? 2.0 : 2.0 / (1.0 + r)) :
            (r < 2.0 / (1.0 + r) ? r : 2.0 / (1.0 + r));
    }
}

#endif
