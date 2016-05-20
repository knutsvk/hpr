#ifndef __HPR_CPP
#define __HPR_CPP

#include "HPR.h"

/* BASE CLASS FOR HPR MATERIALS */

void HyperbolicPeshkovRomenski::xFlux( 
        const SimpleArray< double, 14 >& Q, 
        SimpleArray< double, 14 >& F )
{
    double rho = getDensity( Q );
    SimpleArray< double, 3 > u = getVelocity( Q );
    Eigen::Matrix3d A = getDistortion( Q );
    double E = getEnergy( Q );

    Eigen::Matrix3d sigma = getShearStress( Q );
    double p = getPressure( Q );
    
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

void HyperbolicPeshkovRomenski::yFlux( 
        const SimpleArray< double, 14 >& Q, 
        SimpleArray< double, 14 >& G )
{
    double rho = getDensity( Q );
    SimpleArray< double, 3 > u = getVelocity( Q );
    Eigen::Matrix3d A = getDistortion( Q );
    double E = getEnergy( Q );

    Eigen::Matrix3d sigma = getShearStress( Q );
    double p = getPressure( Q );
    
    G[0] = rho * u[1];
    G[1] = rho * u[1] * u[0] - sigma(0, 1) + p;
    G[2] = rho * u[1] * u[1] - sigma(1, 1);
    G[3] = rho * u[1] * u[2] - sigma(2, 1);
    G[4] = 0.0;
    G[5] = A( 0, 0 ) * u[0] + A(0, 1) * u[1] + A(0, 2) * u [2];
    G[6] = 0.0;
    G[7] = 0.0;
    G[8] = A(1, 0) * u[0] + A(1, 1) * u[1] + A(1, 2) * u [2];
    G[9] = 0.0;
    G[10] = 0.0;
    G[11] = A(2, 0) * u[0] + A(2, 1) * u[1] + A(2, 2) * u [2];
    G[12] = 0.0;
    G[13] = rho * u[1] * E + u[1] * p - sigma(0, 1) * u[0] 
        - sigma(1, 1) * u[1] - sigma(2, 1) * u[2];
}

HyperbolicPeshkovRomenski::HyperbolicPeshkovRomenski( 
        double _shearSoundSpeed, double _referenceDensity, 
        int _nCellsX, int nCellsY, double _domain[4] )
{
    c_s = _shearSoundSpeed; 
    rho_0 = _referenceDensity;
    nCellsX = _nCellsX;
    nCellsY = _nCellsY;

    for(unsigned i = 0; i < 4; i++)
    {
        domain[i] = _domain[i];
    }

    dx = ( domain[1] - domain[0] ) / nCellsX;
    dy = ( domain[3] - domain[2] ) / nCellsY; 

    consVars.resize( ( nCellsX + 2 * nGhostCells ) * ( nCellsY + 2 * nGhostCells ) );
    xDirFlux.resize( ( nCellsX + 1 ) * nCellsY );
    yDirFlux.resize( nCellsX * ( nCellsY + 1 ) );
}

double HyperbolicPeshkovRomenski::getDensity( 
        const SimpleArray< double, 14 >& Q )
{
    return Q[0];
}

SimpleArray< double, 3 > HyperbolicPeshkovRomenski::getVelocity( 
        const SimpleArray< double, 14 >& Q )
{
    SimpleArray< double, 3 > u; 
    for( unsigned i = 0; i < 3; i++ ) 
    { 
        u[i] = Q[i + 1];
    }
    return u / Q[0];
}

Eigen::Matrix3d HyperbolicPeshkovRomenski::getDistortion( 
        const SimpleArray< double, 14 >& Q )
{
    Eigen::Matrix3d A;
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

double HyperbolicPeshkovRomenski::mesoEnergy( Eigen::Matrix3d distortion )
{
    Eigen::Matrix3d G = distortion.trace() * distortion;
    return 0.25 * c_s * c_s * ( G.trace() - ( G * G ).trace() / 3.0 );
}

double HyperbolicPeshkovRomenski::macroEnergy( SimpleArray< double, 3 > u )
{
    return 0.5 * ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
}

void HyperbolicPeshkovRomenski::initialize( double initDiscontRad, 
        double density_L, double density_R, 
        SimpleArray< double, 3 > velocity_L, SimpleArray< double, 3 > velocity_R, 
        Eigen::Matrix3d distortion_L, Eigen::Matrix3d distortion_R,
        double pressure_L, double pressure_R )
{
    double x, y, r;
    int cell;

#pragma omp parallel for private( x, y, r, cell )
    for( unsigned i = 0; i < nCellsX + 2 * nGhostCells; i++ )
    {
        for( unsigned j = 0; j < nCellsY + 2 * nGhostCells; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            x = domain[0] + ( i - nGhostCells + 0.5 ) * dx;
            y = domain[2] + ( j - nGhostCells + 0.5 ) * dy; 
            r = sqrt( x * x + y * y );

            if(r < initDiscontRad)
            {
                consVars[cell][0] = density_L;
                for( unsigned k = 0; k < 3; k++ )
                {
                    consVars[cell][k + 1] = density_L * velocity_L[k];
                    for( unsigned l = 0; l < 3; l++ )
                    {
                        consVars[cell][4 + 3 * k + l] = distortion_L(k, l);
                    }
                }
                consVars[cell][13] = density_L * 
                    ( microEnergy( density_L, pressure_L ) +
                      mesoEnergy( distortion_L ) + 
                      macroEnergy( velocity_L ) );
            }
            else
            {
                consVars[cell][0] = density_R;
                for( unsigned k = 0; k < 3; k++ )
                {
                    consVars[cell][k + 1] = density_R * velocity_R[k];
                    for( unsigned l = 0; l < 3; l++ )
                    {
                        consVars[cell][4 + 3 * k + l] = distortion_R(k, l);
                    }
                }
                consVars[cell][13] = density_R * 
                    ( microEnergy( density_R, pressure_R ) +
                      mesoEnergy( distortion_R ) + 
                      macroEnergy( velocity_R ) );
            }
        }
    }
}

void HyperbolicPeshkovRomenski::transmissiveBCs()
{
    // TODO: ALlow different BCs at each edge
    int nCells = ( nCellsX + 2 * nGhostCells ) * ( nCellsY + 2 * nGhostCells );
    int cell;

    for( unsigned i = 0; i < nGhostCells; i++ )
    {
        for( unsigned j = nGhostCells; j < nGhostCells + nCellsY )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            consVars[cell] = consVars[cell + nCellsY + 2 * nGhostCells];
            consVars[nCells - 1 - cell] 
                = consVars[nCells - 1 - ( cell + nCellsY + 2 * nGhostCells )];
        }
    }

    for( unsigned i = nGhostCells; i < nGhostCells + nCellsX ; i++ )
    {
        for( unsigned j = 0; j < nGhostCells; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            consVars[cell] = consVars[cell + 1];
            consVars[nCells - 1 - cell] = consVars[nCells - 1 - ( cell + 1 )];
        }
    }
}

void HyperbolicPeshkovRomenski::reflectiveBCs()
{
    // TODO
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
    double rho;
    Eigen::Matrix3d A; 
    double scaleFactor; 
#pragma omp parallel for private( rho, A, scaleFactor )
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        rho = getDensity( consVars[i] );
        A = getDistortion( consVars[i] );
        scaleFactor = pow( rho / ( rho_0 * A.determinant() ), 1.0 / 3.0 );
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
    double rho; 
    SimpleArray< double, 3 > u;
    double p; 
    Eigen::Matrix3d sigma;

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

/* CLASS HPR_FLUID */

HPR_Fluid::HPR_Fluid( double _shearSoundSpeed, double _referenceDensity, 
        int _nCells, double _domain[2], 
        double _gamma, double _strainDissipationTime ) : 
    HyperbolicPeshkovRomenski( _shearSoundSpeed, _referenceDensity, 
            _nCells, _domain) 
{
    gamma = _gamma;
    tau = _strainDissipationTime;
}

double HPR_Fluid::getPressure( const SimpleArray< double, 14 >& Q )
{
    double rho = getDensity( Q );
    SimpleArray< double, 3 > u = getVelocity( Q );
    Eigen::Matrix3d A = getDistortion( Q );
    double E = getEnergy( Q );
    
    double E_2 = mesoEnergy( A );
    double E_3 = macroEnergy( u );

    return ( gamma - 1.0 ) * rho * ( E - E_2 - E_3 );
}

Eigen::Matrix3d HPR_Fluid::getSource( const SimpleArray< double, 14 >& Q )
{
    double nu = getDensity( Q ) / rho_0;
    Eigen::Matrix3d A = getDistortion( Q );
    Eigen::Matrix3d G = A.transpose() * A;
    Eigen::Matrix3d devG = G - G.trace() * Eigen::Matrix3d::Identity() / 3.0;
    return - 3.0 * pow(nu, 5.0 / 3.0) * A * devG / tau;
}

double HPR_Fluid::microEnergy( double density, double pressure )
{
    return pressure / ( ( gamma - 1.0 ) * density );
}

double HPR_Fluid::getTimeStep( const double c_CFL )
{
    std::vector< double > S( nCells + 2 * nGhostCells ); 
    double rho;
    SimpleArray< double, 3 > u; 
    double p;

#pragma omp parallel for private( rho, u, p )
    // TODO: reduce for loop with max
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        rho = getDensity( consVars[i] );
        u = getVelocity( consVars[i] );
        p = getPressure( consVars[i] );
        S[i] = fabs( u[0] ) + sqrt( gamma * p / rho + 4.0 * c_s * c_s / 3.0 );
    }

    return c_CFL * dx / *std::max_element( S.begin(), S.end() );
}

void HPR_Fluid::integrateODE( double dt )
{

#pragma omp parallel for 
    for( unsigned i = nGhostCells; i < nCells + nGhostCells; i++ )
    {
        integrate_adaptive( make_controlled( 1.0e-12, 1.0e-12, stepper_type() ), 
                HPR_ODE(tau, rho_0), consVars[i], 0.0, 0.0 + dt, 1.0e-3 * dt );
    }
}

/* CLASS HPR_SOLID */

HPR_Solid::HPR_Solid( double _shearSoundSpeed, double _referenceDensity, 
        int _nCells, double _domain[2], 
        double _c_0, double _Gamma_0, double _s_H ) : 
    HyperbolicPeshkovRomenski( _shearSoundSpeed, _referenceDensity, 
            _nCells, _domain) 
{
    c_0 = _c_0;
    Gamma_0 = _Gamma_0;
    s_H = _s_H;
}

double HPR_Solid::getPressure( const SimpleArray< double, 14 >& Q )
{
    double rho = getDensity( Q );
    SimpleArray< double, 3 > u = getVelocity( Q );
    Eigen::Matrix3d A = getDistortion( Q );
    double E = getEnergy( Q );
    
    double E_2 = mesoEnergy( A );
    double E_3 = macroEnergy( u );

    double nu = rho / rho_0;
    double f = ( nu - 1.0 ) * ( nu - 0.5 * Gamma_0 * ( nu - 1.0 ) ) / 
        pow( ( nu - s_H * ( nu - 1.0 ) ), 2.0 );

    return rho_0 * c_0 * c_0 * f + rho_0 * Gamma_0 * ( E  - E_2 - E_3 );
}

double HPR_Solid::microEnergy( double density, double pressure )
{
    double nu = density / rho_0; 
    double f = ( nu - 1.0 ) * ( nu - 0.5 * Gamma_0 * ( nu - 1.0 ) ) / 
        pow( ( nu - s_H * ( nu - 1.0 ) ), 2.0 );

    return ( pressure - rho_0 * c_0 * c_0 * f ) / ( rho_0 * Gamma_0 );
}

double HPR_Solid::getTimeStep( const double c_CFL )
{
    std::vector< double > S( nCells + 2 * nGhostCells ); 
    SimpleArray< double, 3 > u; 

#pragma omp parallel for private( u )
    // TODO: reduce for loop with max
    for( unsigned i = 0; i < nCells + 2 * nGhostCells; i++ )
    {
        u = getVelocity( consVars[i] );
        S[i] = fabs( u[0] ) + sqrt( c_0 * c_0 + 4.0 * c_s * c_s / 3.0 );
    }

    return c_CFL * dx / *std::max_element( S.begin(), S.end() );
}

/* ODE STRUCT */

HPR_ODE::HPR_ODE( double _tau, double _rho_0 )
{
    tau = _tau; 
    rho_0 = _rho_0;
}

void HPR_ODE::operator()( const state_type& Q, state_type& S, const double t )
{
    double nu = Q[0] / rho_0;

    Eigen::Matrix3d A;
    for( unsigned i = 0; i < 3; i++ )
    {
        for( unsigned j = 0; j < 3; j++ )
        {
            A(i, j) = Q[4 + 3 * i + j];
        }
    }

    Eigen::Matrix3d G = A.transpose() * A;
    Eigen::Matrix3d devG = G - G.trace() / 3.0 * Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Psi = - 3.0 * pow( nu, 5.0 / 3.0 ) / tau * A * devG;

    // TODO: std::fill to set elements of S to zero by default
    for( unsigned i = 0; i < 4; i++ )
    { // no source terms for density, velocity
        S[i] = 0.0;
    }

    for( unsigned i = 0; i < 3; i++ )
    {
        for( unsigned j = 0; j < 3; j++ )
        { // source terms for distortion tensor
            S[4 + 3 * i + j] = Psi(i, j);
        }
    }

    S[13] = 0.0; // no source term for energy
}

/* CLASS INDPDNT FUNCTIONS */

double slopeLimiter( double q_min, double q_0, double q_plus )
{ // TODO: move to different file Limiters.cpp, allow minbee vanleer etc
    double Delta_L, Delta_R; 

    if( fabs( q_0 - q_min ) < 1.0e-16 ) 
        Delta_L = copysign( 1.0e-16, q_0 - q_min );
    else 
        Delta_L = q_0 - q_min;

    if( fabs( q_plus -q_0 ) < 1.0e-16 )
        Delta_R = copysign( 1.0e-16, q_plus - q_0 );
    else 
        Delta_R = q_plus -q_0;
    
    double r = Delta_L / Delta_R;

    // superbee
    if( r <= 0.0 )
        return 0.0;
    else if( r <= 0.5 ) 
        return 2.0 * r; 
    else if( r <= 1.0 ) 
        return 1.0; 
    else //TODO: std::min_element()
        return ( 2.0 < r ) ?  ( 2.0 < 2.0 / ( 1.0 + r ) ? 2.0 : 2.0 / ( 1.0 + r ) ) :
            ( r < 2.0 / ( 1.0 + r ) ? r : 2.0 / ( 1.0 + r ) );
}

#endif
