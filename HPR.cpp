#ifndef __HPR_CPP
#define __HPR_CPP

#include "HPR.h"

/* BASE CLASS FOR EULER SOLVER */

void Euler::xFlux(
        const SimpleArray< double, 4 >& Q,
        SimpleArray< double, 4 >& F )
{
    double rho = getDensity( Q );
    SimpleArray< double, 2 > u = getVelocity( Q );
    double E = getEnergy( Q );
    double p = getPressure( Q );

    F[0] = rho * u[0];
    F[1] = rho * u[0] * u[0] + p;
    F[2] = rho * u[0] * u[1];
    F[3] = u[0] * ( rho * E + p );
}

void Euler::yFlux(
        const SimpleArray< double, 4 >& Q,
        SimpleArray< double, 4 >& G )
{
    double rho = getDensity( Q );
    SimpleArray< double, 2 > u = getVelocity( Q );
    double E = getEnergy( Q );
    double p = getPressure( Q );

    G[0] = rho * u[1];
    G[1] = rho * u[1] * u[0];
    G[2] = rho * u[1] * u[1] + p;
    G[3] = u[1] * ( rho * E + p );
}

void Euler::forceFlux( double dt, double dr, int dir,
        const SimpleArray< double, 4 >& Q_L,
        const SimpleArray< double, 4 >& Q_R,
        SimpleArray< double, 4 >& F )
{
    SimpleArray< double, 4 > F_L, F_R, Q_0, F_0;

    if( dir == 0 )
    {
        xFlux( Q_L, F_L );
        xFlux( Q_R, F_R );
    }
    else
    {
        yFlux( Q_L, F_L );
        yFlux( Q_R, F_R );
    }

    Q_0 = 0.5 * ( Q_L + Q_R ) + 0.5 * dt / dr * ( F_L - F_R );

    if( dir == 0 )
        xFlux( Q_0, F_0 );
    else
        yFlux( Q_0, F_0 );

    F = 0.5 * ( F_0 + 0.5 * ( F_L + F_R ) ) + 0.25 * dr / dt * ( Q_L - Q_R );
}

void Euler::slicFlux( double dt, double dr, int dir,
        const SimpleArray< double, 4 >& Q_2L,
        const SimpleArray< double, 4 >& Q_L,
        const SimpleArray< double, 4 >& Q_R,
        const SimpleArray< double, 4 >& Q_2R,
        SimpleArray< double, 4 >& F )
{
    SimpleArray< double, 4 > xi_L, xi_R,
        Q_L_plus, Q_R_plus, Q_L_0, Q_R_0,
        F_L_plus, F_R_plus, F_L_0, F_R_0,
        Q_L_bar, Q_R_bar;

    // Calculate TVD slope limiters
    for( int j = 0; j < 4; j++ )
    {
        xi_L[j] = slopeLimiter( Q_2L[j], Q_L[j], Q_R[j] );
        xi_R[j] = slopeLimiter( Q_L[j], Q_R[j], Q_2R[j] );
    }

    // Boundary extrapolated values
    Q_L_plus = Q_R - 0.25 * xi_R * ( Q_2R - Q_L );
    Q_R_plus = Q_R + 0.25 * xi_R * ( Q_2R - Q_L );
    Q_L_0 = Q_L - 0.25 * xi_L * ( Q_R - Q_2L );
    Q_R_0 = Q_L + 0.25 * xi_L * ( Q_R - Q_2L );

    if( dir == 0 )
    {
        xFlux( Q_L_plus, F_L_plus );
        xFlux( Q_R_plus, F_R_plus );
        xFlux( Q_L_0, F_L_0 );
        xFlux( Q_R_0, F_R_0 );
    }
    else if( dir == 1 )
    {
        yFlux( Q_L_plus, F_L_plus );
        yFlux( Q_R_plus, F_R_plus );
        yFlux( Q_L_0, F_L_0 );
        yFlux( Q_R_0, F_R_0 );
    }

    // Evolve by time 0.5 * dt
    Q_L_bar = Q_L_plus + 0.5 * dt / dr * ( F_L_plus - F_R_plus );
    Q_R_bar = Q_R_0 + 0.5 * dt / dr * ( F_L_0 - F_R_0 );

    forceFlux( dt, dr, dir, Q_R_bar, Q_L_bar, F );
}

Euler::Euler( int _nCellsX, int _nCellsY, double _domain[4], double _gamma )
{
    nCellsX = _nCellsX;
    nCellsY = _nCellsY;
    nGhostCells = 2;
    nCellsTot = ( nCellsX + 2 * nGhostCells ) * ( nCellsY + 2 * nGhostCells );

    for(int i = 0; i < 4; i++)
        domain[i] = _domain[i];

    dx = ( domain[1] - domain[0] ) / nCellsX;
    dy = ( domain[3] - domain[2] ) / nCellsY;

    gamma = _gamma;

    consVars.resize( nCellsTot );
}

double Euler::getDensity( const SimpleArray< double, 4 >& Q )
{
    return Q[0];
}

SimpleArray< double, 2 > Euler::getVelocity( const SimpleArray< double, 4 >& Q)
{
    SimpleArray< double, 2 > u;
    for( int i = 0; i < 2; i++ )
        u[i] = Q[i + 1] / Q[0];

    return u;
}

double Euler::getEnergy(
        const SimpleArray< double, 4 >& Q )
{
    return Q[3] / Q[0];
}

double Euler::getPressure( const SimpleArray< double, 4 >& Q )
{
    double rho = getDensity( Q );
    SimpleArray< double, 2 > u = getVelocity( Q );
    double E = getEnergy( Q );
    double E_mac = macroEnergy( u );

    return ( gamma - 1.0 ) * rho * ( E - E_mac );
}

double Euler::microEnergy( double density, double pressure )
{
    return pressure / ( ( gamma - 1.0 ) * density );
}

double Euler::macroEnergy( SimpleArray< double, 2 > u )
{
    return 0.5 * ( u[0] * u[0] + u[1] * u[1] );
}

double Euler::getTimeStep( const double c_CFL )
{
    std::vector< double > Sx( nCellsTot );
    std::vector< double > Sy( nCellsTot );
    double rho;
    SimpleArray< double, 2 > u;
    double p;
    double a;

#pragma omp parallel for private( rho, u, p, a )
    // TODO: reduce for loop with max
    for( int i = 0; i < nCellsTot; i++ )
    {
        rho = getDensity( consVars[i] );
        u = getVelocity( consVars[i] );
        p = getPressure( consVars[i] );
        a = sqrt( gamma * p / rho );
        Sx[i] = fabs( u[0] ) + a;
        Sy[i] = fabs( u[1] ) + a;
    }
    double SxMax = *std::max_element( Sx.begin(), Sx.end() );
    double SyMax = *std::max_element( Sy.begin(), Sy.end() );
    return c_CFL * std::min( dx / SxMax, dy / SyMax);
}

void Euler::initialise( double initDiscontPos, Direction initDiscontDir, 
        double density[2], SimpleArray< double, 2 > velocity[2], 
        double pressure[2] )
{
    int cell;
    double x, y, r;
    bool isLeft = true;

    for( int i = 0; i < nCellsX + 2 * nGhostCells; i++ )
    {
        for( int j = 0; j < nCellsY + 2 * nGhostCells; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            x = domain[0] + ( i - nGhostCells + 0.5 ) * dx;
            y = domain[2] + ( j - nGhostCells + 0.5 ) * dy;
            r = sqrt( x * x + y * y );

            switch( initDiscontDir )
            {
                case horizontal:
                    isLeft = x <= initDiscontPos;
                    break;
                case vertical:
                    isLeft = y <= initDiscontPos;
                    break;
                case radial:
                    isLeft = r <= initDiscontPos;
                    break;
            }

            if( isLeft )
            {
                consVars[cell][0] = density[0];
                for( int k = 0; k < 2; k++ )
                {
                    consVars[cell][k + 1] = density[0] * velocity[0][k];
                }
                consVars[cell][3] = density[0] * 
                    ( microEnergy( density[0], pressure[0] ) + 
                      macroEnergy( velocity[0] ) );
            }
            else
            {
                consVars[cell][0] = density[1];
                for( int k = 0; k < 2; k++ )
                {
                    consVars[cell][k + 1] = density[1] * velocity[1][k];
                }
                consVars[cell][3] = density[1] *
                    ( microEnergy( density[1], pressure[1] ) +
                      macroEnergy( velocity[1] ) );
            }
        }
    }
}

void Euler::initialiseConvergenceTest()
{
    int cell;
    const double EPS = 5.0;
    const double GAM = 1.4;
    double x, y, r2, rho, p;
    SimpleArray< double, 2 > u;

    for( int i = 0; i < nCellsX + 2 * nGhostCells; i++ )
    {
        for( int j = 0; j < nCellsY + 2 * nGhostCells; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            x = domain[0] + ( i - nGhostCells + 0.5 ) * dx;
            y = domain[2] + ( j - nGhostCells + 0.5 ) * dy;
            r2 = x * x + y * y;
            rho = pow( 1.0 - ( GAM - 1.0 ) * EPS * EPS / ( 8 * GAM * M_PI *
                        M_PI ) * exp( 1 - r2 ), 1.0 / ( GAM - 1.0 ) );
            u[0] = 1.0 - y * EPS / ( 2 * M_PI) * exp( 0.5 * ( 1.0 - r2 ) );
            u[1] = 1.0 + x * EPS / ( 2 * M_PI) * exp( 0.5 * ( 1.0 - r2 ) );
            p = pow( rho, GAM );

            consVars[cell][0] = rho;
            consVars[cell][1] = rho * u[0];
            consVars[cell][2] = rho * u[1];
            consVars[cell][3] = rho * ( microEnergy( rho, p ) + 
                    macroEnergy( u) );
        }
    }
}

void Euler::exactConvergenceSolution()
{
    int cell;
    char outfile[100];
    const double EPS = 5.0;
    const double GAM = 1.4;
    const double T_F = 1.0;
    double x, y, r2, rho, p;
    SimpleArray< double, 2 > u;

    for( int i = 0; i < nCellsX + 2 * nGhostCells; i++ )
    {
        for( int j = 0; j < nCellsY + 2 * nGhostCells; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            x = domain[0] + ( i - nGhostCells + 0.5 ) * dx - T_F;
            y = domain[2] + ( j - nGhostCells + 0.5 ) * dy - T_F;
            r2 = x * x + y * y;
            rho = pow( 1.0 - ( GAM - 1.0 ) * EPS * EPS / ( 8 * GAM * M_PI *
                        M_PI ) * exp( 1 - r2 ), 1.0 / ( GAM - 1.0 ) );
            u[0] = 1.0 - y * EPS / ( 2 * M_PI) * exp( 0.5 * ( 1.0 - r2 ) );
            u[1] = 1.0 + x * EPS / ( 2 * M_PI) * exp( 0.5 * ( 1.0 - r2 ) );
            p = pow( rho, GAM );

            consVars[cell][0] = rho;
            consVars[cell][1] = rho * u[0];
            consVars[cell][2] = rho * u[1];
            consVars[cell][3] = rho * ( microEnergy( rho, p )
                + macroEnergy( u ) );
        }
    }
    sprintf( outfile, "./Results/ConvergenceStudies_1DX_Nx%d_Ny%d_Exact.out",
            nCellsX, nCellsY );
    output1DSliceX( outfile );
    sprintf( outfile, "./Results/ConvergenceStudies_1DY_Nx%d_Ny%d_Exact.out",
            nCellsX, nCellsY );
    output1DSliceY( outfile );
    sprintf( outfile, "./Results/ConvergenceStudies_2D_Nx%d_Ny%d_Exact.out",
            nCellsX, nCellsY );
    output2D( outfile );
}

void Euler::boundaryConditions( BoundaryCondition type[4] )
{
    /* type[0]: left BC type
     * type[1]: right BC type
     * type[2]: bottom BC type
     * type[3]: top BC type
     */

    int cell, copyTo, copyFrom;
    int M = nCellsY + 2 * nGhostCells;

    for( int i = 0; i < nGhostCells; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * M + j;

            // Left boundary (x = xMin)
            switch ( type[0] )
            {
                case transmissive:
                    copyTo = cell;
                    copyFrom = copyTo + M;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyTo = cell;
                    copyFrom = copyTo + ( 2 * ( nGhostCells - i ) - 1 ) * M;
                    consVars[copyTo] = consVars[copyFrom];
                    for( int k = 1; k < 3; k++ )
                    {
                        consVars[copyTo][k] *= - 1.0;
                    }
                    break;
                case periodic:
                    copyTo = cell;
                    copyFrom = cell + nCellsX * M;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }

            // Right boundary (x = xMax)
            switch (type[1] )
            {
                case transmissive:
                    copyTo = nCellsTot - ( 1 + cell );
                    copyFrom = copyTo - M;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyTo = nCellsTot - ( 1 + cell );
                    copyFrom = copyTo - ( 2 * ( nGhostCells - i ) - 1 ) * M;
                    consVars[copyTo] = consVars[copyFrom];
                    for( int k = 1; k < 3; k++ )
                    {
                        consVars[copyTo][k] *= - 1.0;
                    }
                    break;
                case periodic:
                    copyTo = cell + ( nCellsX + nGhostCells ) * M;
                    copyFrom = cell + nGhostCells * M;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }
        }
    }

    for( int i = nGhostCells; i < nGhostCells + nCellsX ; i++ )
    {
        for( int j = 0; j < nGhostCells; j++ )
        {
            cell = i * M + j;

            // Bottom boundary (y = yMin)
            switch( type[2] )
            {
                case transmissive:
                    copyTo = cell;
                    copyFrom = copyTo + 1;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyTo = cell;
                    copyFrom = copyTo + ( 2 * ( nGhostCells - j ) - 1 );
                    consVars[copyTo] = consVars[copyFrom];
                    for( int k = 1; k < 3; k++ )
                    {
                        consVars[copyTo][k] *= - 1.0;
                    }
                    break;
                case periodic:
                    copyTo = cell;
                    copyFrom = cell + nCellsY;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }

            // Top boundary (y = yMax)

            switch( type[3] )
            {
                case transmissive:
                    copyTo = ( i + 1 ) * M - ( 1 + j );
                    copyFrom = copyTo - 1;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyTo = ( i + 1 ) * M - ( 1 + j );
                    copyFrom = copyTo - ( 2 * ( nGhostCells - j ) - 1 );
                    consVars[copyTo] = consVars[copyFrom];
                    for( int k = 1; k < 3; k++ )
                    {
                        consVars[copyTo][k] *= - 1.0;
                    }
                    consVars[copyTo][1] += 2.0;
                    break;
                case periodic:
                    copyTo = cell + nCellsY + nGhostCells;
                    copyFrom = cell + nGhostCells;
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }
        }
    }

    if( !isPhysical() )
    {
        std::cout << "Unphysical state encountered in function "
            << "boundaryConditions() " << std::endl;
    }
}

void Euler::xSweep( double dt )
{
    std::vector< SimpleArray< double, 4 > > tempVars( nCellsTot );
#pragma omp parallel for
    for( int i = 0; i < nCellsTot; i++ )
        tempVars[i] = consVars[i];

    int cell;
    int M = nCellsY + 2 * nGhostCells;
    SimpleArray< double, 4 > F_L, F_R;

#pragma omp parallel for private( cell, F_L, F_R )
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * M + j;
            slicFlux( dt, dx, 0,
                    tempVars[cell - 2 * M],
                    tempVars[cell - M],
                    tempVars[cell],
                    tempVars[cell + M],
                    F_L );
            slicFlux( dt, dx, 0,
                    tempVars[cell - M],
                    tempVars[cell],
                    tempVars[cell + M],
                    tempVars[cell + 2 * M],
                    F_R );
            consVars[cell] = tempVars[cell] + dt / dx * ( F_L - F_R );
        }
    }

    if( !isPhysical() )
    {
        std::cout << "Unphysical state encountered in function "
            << "xSweep() " << std::endl;
    }
}

void Euler::ySweep( double dt )
{
    std::vector< SimpleArray< double, 4 > > tempVars( nCellsTot );
#pragma omp parallel for
    for( int i = 0; i < nCellsTot; i++ )
        tempVars[i] = consVars[i];

    int cell;
    SimpleArray< double, 4 > F_B, F_T;

#pragma omp parallel for private( cell, F_B, F_T )
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            slicFlux( dt, dy, 1,
                    tempVars[cell - 2],
                    tempVars[cell - 1],
                    tempVars[cell],
                    tempVars[cell + 1],
                    F_B );
            slicFlux( dt, dy, 1,
                    tempVars[cell - 1],
                    tempVars[cell],
                    tempVars[cell + 1],
                    tempVars[cell + 2],
                    F_T );
            consVars[cell] = tempVars[cell] + dt / dy * ( F_B - F_T );
        }
    }

    if( !isPhysical() )
    {
        std::cout << "Unphysical state encountered in function "
            << "ySweep() " << std::endl;
    }
}

void Euler::output2D( char* filename )
{
    int M = nCellsY + 2 * nGhostCells;
    int cell;
    double x, y;
    double rho, p, e;
    SimpleArray< double, 2 > u;

    std::ofstream fs;
    fs.open( filename );
    fs << "x" << "\t" << "y" << "\t" << "rho" << "\t" << "u" << "\t" << "v" 
        << "\t" << "p" << "\t" << "e" << std::endl;

    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * M + j;
            x = domain[0] + (i - nGhostCells + 0.5 ) * dx;
            y = domain[2] + (j - nGhostCells + 0.5 ) * dy;
            rho = getDensity( consVars[cell] );
            u = getVelocity( consVars[cell] );
            p = getPressure( consVars[cell] );
            e = microEnergy( rho, p) ;

            fs << x << " \t" << y << " \t" << rho << " \t" << u[0] << " \t" 
                << u[1] << " \t" << p << "\t" << e << std::endl;
        }
        fs << std::endl;
    }
    fs.close();
}

void Euler::output1DSliceX( char* filename )
{
    int cell;
    double x;
    double rho, p, e;
    SimpleArray< double, 2 > u;

    std::ofstream fs;
    fs.open( filename );
    fs << "x" << "\t" << "rho" << "\t" << "u" << "\t" << "v" << "\t" << "p" 
        << "\t" << "e" << std::endl;

    int j = nGhostCells + nCellsY / 2;
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        cell = i * ( nCellsY + 2 * nGhostCells ) + j;
        x = domain[0] + (i - nGhostCells + 0.5 ) * dx;
        rho = getDensity( consVars[cell] );
        u = getVelocity( consVars[cell] );
        p = getPressure( consVars[cell] );
        e = microEnergy( rho, p );

        fs << x << "\t" << rho << "\t" << u[0] << "\t" << u[1] << "\t"
            << p << "\t" << e << std::endl;
    }
    fs.close();
}

void Euler::output1DSliceY( char* filename )
{
    int cell;
    double y;
    double rho, p, e;
    SimpleArray< double, 2 > u;

    std::ofstream fs;
    fs.open( filename );
    fs << "y" << "\t" << "rho" << "\t" << "u" << "\t" << "v" << "\t"
        << "p" << "\t" << "e" << std::endl;

    int i = nGhostCells + nCellsX / 2;
    for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
    {
        cell = i * ( nCellsY + 2 * nGhostCells ) + j;
        y = domain[2] + (j - nGhostCells + 0.5 ) * dy;
        rho = getDensity( consVars[cell] );
        u = getVelocity( consVars[cell] );
        p = getPressure( consVars[cell] );
        e = microEnergy( rho, p );

        fs << y << "\t" << rho << "\t" << u[0] << "\t" << u[1] << "\t"
            << p << "\t" << e << std::endl;
    }
    fs.close();
}

bool Euler::isPhysical()
{
    int M = nCellsY + 2 * nGhostCells;
    double x, y, rho, E, p;
    SimpleArray< double, 2 > u;
    for( int i = 0; i < nCellsX + 2 * nGhostCells; i++ )
    {
        for( int j = 0; j < nCellsY + 2 * nGhostCells; j++ )
        {
            x = domain[0] + ( i - nGhostCells + 0.5 ) * dx;
            y = domain[2] + ( j - nGhostCells + 0.5 ) * dy;

            rho = getDensity( consVars[i * M + j] );
            if( std::isnan( rho ) )
            {
                std::cout << "Error: rho = nan at x = " << x << ", y = " << y
                    << ". " << std::endl;
                return false;
            }
            else if( rho < 0.0 )
            {
                std::cout << "Error: rho < 0.0 at x = " << x << ", y = " << y
                    << ". " << std::endl;
                return false;
            }

            u = getVelocity( consVars[i * M + j] );
            for( int k = 0; k < 2; k++ )
            {
                if( std::isnan( u[k] ) )
                {
                    std::cout << "Error: u[" << k << "] = nan at x = " << x
                        << ", y = " << y << ". " << std::endl;
                    return false;
                }
            }

            E = getEnergy( consVars[i * M + j] );
            if( std::isnan( E ) )
            {
                std::cout << "Error: E = nan at x = " << x << ", y = " << y
                    << ". " << std::endl;
                return false;
            }

            p = getPressure( consVars[i * M + j] );
            if( std::isnan( p ) )
            {
                std::cout << "Error: p = nan at x = " << x << ", y = " << y
                    << ". " << std::endl;
                return false;
            }
        }
    }
    return true;
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

void configurate( const char* inputFile, int& nCellsX, int& nCellsY, 
        double& CFL, double& tStop, double& gamma, double domain[4], 
        double& initDiscontPos, Direction& initDiscontDir, 
        BoundaryCondition BCs[4], double initDensity[2], 
        SimpleArray< double, 2 > initVelocity[2], double initPressure[2] )
{
    // Read input file
    libconfig::Config cfg;
    cfg.readFile( inputFile );

    int BCs_int[4], dir_int;

    cfg.lookupValue( "nCellsX", nCellsX );
    cfg.lookupValue( "nCellsY", nCellsY );
    cfg.lookupValue( "CFL", CFL );
    cfg.lookupValue( "tStop", tStop );
    cfg.lookupValue( "heatCapacityRatio", gamma );
    cfg.lookupValue( "xMin", domain[0] );
    cfg.lookupValue( "xMax", domain[1] );
    cfg.lookupValue( "yMin", domain[2] );
    cfg.lookupValue( "yMax", domain[3] );
    cfg.lookupValue( "initDiscontDir", dir_int );
    cfg.lookupValue( "initDiscontPos", initDiscontPos );
    cfg.lookupValue( "leftBC", BCs_int[0] );
    cfg.lookupValue( "rightBC", BCs_int[1] );
    cfg.lookupValue( "bottomBC", BCs_int[2] );
    cfg.lookupValue( "topBC", BCs_int[3] );
    cfg.lookupValue( "rho_1", initDensity[0] );
    cfg.lookupValue( "u_1", initVelocity[0][0] );
    cfg.lookupValue( "v_1", initVelocity[0][1] );
    cfg.lookupValue( "p_1", initPressure[0] );
    cfg.lookupValue( "rho_2", initDensity[1] );
    cfg.lookupValue( "u_2", initVelocity[1][0] );
    cfg.lookupValue( "v_2", initVelocity[1][1] );
    cfg.lookupValue( "p_2", initPressure[1] );

    for( int i = 0; i < 4; i++ )
    {
        BCs[i] = (BoundaryCondition) BCs_int[i];
    }

    initDiscontDir = (Direction) dir_int;
}

#endif
