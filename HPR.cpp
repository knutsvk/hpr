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
    G[1] = rho * u[1] * u[0] - sigma(0, 1);
    G[2] = rho * u[1] * u[1] - sigma(1, 1) + p;
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

void HyperbolicPeshkovRomenski::forceFlux( double dt, double dr, int dir, 
        const SimpleArray< double, 14 >& Q_L, 
        const SimpleArray< double, 14 >& Q_R, 
        SimpleArray< double, 14 >& F )
{
    SimpleArray< double, 14 > F_L, F_R, Q_0, F_0;

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

    Q_0 = 0.5 * ( Q_L + Q_R )
        + 0.5 * dt / dr * ( F_L - F_R );

    if( dir == 0 )
    {
        xFlux( Q_0, F_0 );
    }
    else 
    {
        yFlux( Q_0, F_0 );
    }

    F = 0.5 * ( F_0 + 0.5 * ( F_L + F_R ) ) 
        + 0.25 * dr / dt * ( Q_L - Q_R ); 
}

void HyperbolicPeshkovRomenski::slicFlux( double dt, double dr, int dir, 
        const SimpleArray< double, 14 >& Q_2L, 
        const SimpleArray< double, 14 >& Q_L, 
        const SimpleArray< double, 14 >& Q_R, 
        const SimpleArray< double, 14 >& Q_2R, 
        SimpleArray< double, 14 >& F )
{
    SimpleArray< double, 14 > xi_L, xi_R, 
        Q_L_plus, Q_R_plus, Q_L_0, Q_R_0,
        F_L_plus, F_R_plus, F_L_0, F_R_0,
        Q_L_bar, Q_R_bar; 

    // Calculate TVD slope limiters
    for( int j = 0; j < 14; j++ )
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
    else 
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

void HyperbolicPeshkovRomenski::finiteDiffNoncons( 
        const SimpleArray< double, 14 >& Q_0, 
        const SimpleArray< double, 14 >& Q_L, 
        const SimpleArray< double, 14 >& Q_R, 
        const SimpleArray< double, 14 >& Q_B, 
        const SimpleArray< double, 14 >& Q_T, 
        SimpleArray< double, 14 >& N )
{
    Eigen::Matrix3d A_L = getDistortion( Q_L );
    Eigen::Matrix3d A_R = getDistortion( Q_R );
    Eigen::Matrix3d A_B = getDistortion( Q_B );
    Eigen::Matrix3d A_T = getDistortion( Q_T );
    SimpleArray< double, 3> u = getVelocity( Q_0 );

    for( int i = 0; i < 4; i++ )
    {
        N[i] = 0.0;
    }

    N[4] = u[1] * ( ( A_R(0, 1) - A_L(0, 1) ) / ( 2.0 * dx ) 
            - ( A_T(0, 0) - A_B(0, 0) ) / ( 2.0 * dy ) );

    N[5] = u[0] * ( ( A_T(0, 0) - A_B(0, 0) ) / ( 2.0 * dy ) 
            - ( A_R(0, 1) - A_L(0, 1) ) / ( 2.0 * dx ) );

    N[6] = - u[0] * ( A_R(0, 2) - A_L(0, 2) ) / ( 2.0 * dx )
        - u[1] * ( A_T(0, 2) - A_B(0, 2) ) / ( 2.0 * dy );

    N[7] = u[1] * ( ( A_R(1, 1) - A_L(1, 1) ) / ( 2.0 * dx ) 
            - ( A_T(1, 0) - A_B(1, 0) ) / ( 2.0 * dy ) );

    N[8] = u[0] * ( ( A_T(1, 0) - A_B(1, 0) ) / ( 2.0 * dy ) 
            - ( A_R(1, 1) - A_L(1, 1) ) / ( 2.0 * dx ) );

    N[9] = - u[0] * ( A_R(1, 2) - A_L(1, 2) ) / ( 2.0 * dx )
        - u[1] * ( A_T(1, 2) - A_B(1, 2) ) / ( 2.0 * dy );

    N[10] = u[1] * ( ( A_R(2, 1) - A_L(2, 1) ) / ( 2.0 * dx ) 
            - ( A_T(2, 0) - A_B(2, 0) ) / ( 2.0 * dy ) );

    N[11] = u[0] * ( ( A_T(2, 0) - A_B(2, 0) ) / ( 2.0 * dy ) 
            - ( A_R(2, 1) - A_L(2, 1) ) / ( 2.0 * dx ) );

    N[12] = - u[0] * ( A_R(2, 2) - A_L(2, 2) ) / ( 2.0 * dx )
        - u[1] * ( A_T(2, 2) - A_B(2, 2) ) / ( 2.0 * dy );

    N[13] = 0.0;
}

void HyperbolicPeshkovRomenski::boundaryExtrapolatedNoncons( 
        const SimpleArray< double, 14 >& Q_0, 
        const SimpleArray< double, 14 >& Q_L, 
        const SimpleArray< double, 14 >& Q_R, 
        const SimpleArray< double, 14 >& Q_B, 
        const SimpleArray< double, 14 >& Q_T, 
        SimpleArray< double, 14 >& N )
{
    SimpleArray< double, 14 > xi_x, xi_y, Q_LI, Q_RI, Q_BI, Q_TI; 

    // Calculate TVD slope limiters
    for( int j = 0; j < 14; j++ )
    {
        xi_x[j] = slopeLimiter( Q_L[j], Q_0[j], Q_R[j] );
        xi_y[j] = slopeLimiter( Q_B[j], Q_0[j], Q_T[j] );
    }

    // Boundary extrapolated values
    Q_LI = Q_0 - 0.25 * xi_x * ( Q_R - Q_L );
    Q_RI = Q_0 + 0.25 * xi_x * ( Q_R - Q_L );
    Q_BI = Q_0 - 0.25 * xi_y * ( Q_T - Q_B );
    Q_TI = Q_0 + 0.25 * xi_y * ( Q_T - Q_B );

    Eigen::Matrix3d A_L = getDistortion( Q_LI );
    Eigen::Matrix3d A_R = getDistortion( Q_RI );
    Eigen::Matrix3d A_B = getDistortion( Q_BI );
    Eigen::Matrix3d A_T = getDistortion( Q_TI );

    for( int i = 0; i < 4; i++ )
    {
        N[i] = 0.0;
    }

    N[4] = u[1] * ( ( A_R(0, 1) - A_L(0, 1) ) / dx 
            - ( A_T(0, 0) - A_B(0, 0) ) / dy );

    N[5] = u[0] * ( ( A_T(0, 0) - A_B(0, 0) ) / dy
            - ( A_R(0, 1) - A_L(0, 1) ) / dx );

    N[6] = - u[0] * ( A_R(0, 2) - A_L(0, 2) ) / dx
        - u[1] * ( A_T(0, 2) - A_B(0, 2) ) / dy;

    N[7] = u[1] * ( ( A_R(1, 1) - A_L(1, 1) ) / dx
            - ( A_T(1, 0) - A_B(1, 0) ) / dy );

    N[8] = u[0] * ( ( A_T(1, 0) - A_B(1, 0) ) / dy
            - ( A_R(1, 1) - A_L(1, 1) ) / dx );

    N[9] = - u[0] * ( A_R(1, 2) - A_L(1, 2) ) / dx
        - u[1] * ( A_T(1, 2) - A_B(1, 2) ) / dy;

    N[10] = u[1] * ( ( A_R(2, 1) - A_L(2, 1) ) / dx 
            - ( A_T(2, 0) - A_B(2, 0) ) / dy );

    N[11] = u[0] * ( ( A_T(2, 0) - A_B(2, 0) ) / dy 
            - ( A_R(2, 1) - A_L(2, 1) ) / dx );

    N[12] = - u[0] * ( A_R(2, 2) - A_L(2, 2) ) / dx
        - u[1] * ( A_T(2, 2) - A_B(2, 2) ) / dy;

    N[13] = 0.0;
}

HyperbolicPeshkovRomenski::HyperbolicPeshkovRomenski( 
        double _shearSoundSpeed, double _referenceDensity, 
        int _nCellsX, int _nCellsY, double _domain[4] )
{
    c_s = _shearSoundSpeed; 
    rho_0 = _referenceDensity;
    nCellsX = _nCellsX;
    nCellsY = _nCellsY;
    nGhostCells = 2;
    nCellsTot = ( nCellsX + 2 * nGhostCells ) * ( nCellsY + 2 * nGhostCells );

    for(int i = 0; i < 4; i++)
    {
        domain[i] = _domain[i];
    }

    dx = ( domain[1] - domain[0] ) / nCellsX;
    dy = ( domain[3] - domain[2] ) / nCellsY; 

    consVars.resize( nCellsTot );
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
    for( int i = 0; i < 3; i++ ) 
    { 
        u[i] = Q[i + 1];
    }
    return u / Q[0];
}

Eigen::Matrix3d HyperbolicPeshkovRomenski::getDistortion( 
        const SimpleArray< double, 14 >& Q )
{
    Eigen::Matrix3d A;
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            A(i, j) = Q[4 + 3 * i + j];
        }
    }
    return A;
}

Eigen::Matrix3d HyperbolicPeshkovRomenski::getCurlTerm( 
        const SimpleArray< double, 14 >& Q_0, 
        const SimpleArray< double, 14 >& Q_L, 
        const SimpleArray< double, 14 >& Q_R,
        const SimpleArray< double, 14 >& Q_B, 
        const SimpleArray< double, 14 >& Q_T )
{
    Eigen::Matrix3d A_L = getDistortion( Q_L );
    Eigen::Matrix3d A_R = getDistortion( Q_R );
    Eigen::Matrix3d A_B = getDistortion( Q_B );
    Eigen::Matrix3d A_T = getDistortion( Q_T );
    SimpleArray< double, 3 > u = getVelocity( Q_0 );
    Eigen::Matrix3d B; 

    B(0, 0) = u[1] * ( ( A_T(0, 0) - A_B(0, 0) ) / ( 2.0 * dy ) 
        - ( A_R(0, 1) - A_L(0, 1) ) / ( 2.0 * dx ) );

    B(0, 1) = u[0] * ( ( A_R(0, 1) - A_L(0, 1) ) / ( 2.0 * dx ) 
        - ( A_T(0, 0) - A_B(0, 0) ) / ( 2.0 * dy ) );

    B(0, 2) = u[0] * ( ( A_R(0, 2) - A_L(0, 2) ) / ( 2.0 * dx ) ) 
        + u[1] * ( ( A_T(0, 2) - A_B(0, 2) ) / ( 2.0 * dy ) );

    B(1, 0) = u[1] * ( ( A_T(1, 0) - A_B(1, 0) ) / ( 2.0 * dy ) 
        - ( A_R(1, 1) - A_L(1, 1) ) / ( 2.0 * dx ) );

    B(1, 1) = u[0] * ( ( A_R(1, 1) - A_L(1, 1) ) / ( 2.0 * dx ) 
        - ( A_T(1, 0) - A_B(1, 0) ) / ( 2.0 * dy ) );

    B(1, 2) = u[0] * ( ( A_R(1, 2) - A_L(1, 2) ) / ( 2.0 * dx ) ) 
        + u[1] * ( ( A_T(1, 2) - A_B(1, 2) ) / ( 2.0 * dy ) );

    B(2, 0) = u[1] * ( ( A_T(2, 0) - A_B(2, 0) ) / ( 2.0 * dy ) 
        - ( A_R(2, 1) - A_L(2, 1) ) / ( 2.0 * dx ) );

    B(2, 1) = u[0] * ( ( A_R(2, 1) - A_L(2, 1) ) / ( 2.0 * dx ) 
        - ( A_T(2, 0) - A_B(2, 0) ) / ( 2.0 * dy ) );

    B(2, 2) = u[0] * ( ( A_R(2, 2) - A_L(2, 2) ) / ( 2.0 * dx ) ) 
        + u[1] * ( ( A_T(2, 2) - A_B(2, 2) ) / ( 2.0 * dy ) );

    return B;
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
    return 0.25 * c_s * c_s * ( ( G * G ).trace() - G.trace() * G.trace() / 3.0 );
}

double HyperbolicPeshkovRomenski::macroEnergy( SimpleArray< double, 3 > u )
{
    return 0.5 * ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
}

void HyperbolicPeshkovRomenski::initialize( double initDiscontPos, 
        Direction initDiscontDir, double density[2], 
        SimpleArray< double, 3 > velocity[2], Eigen::Matrix3d distortion[2],
        double pressure[2] )
{
    int cell;
    double x, y, r;
    bool isLeft;

#pragma omp parallel for private( x, y, r, cell )
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
                for( int k = 0; k < 3; k++ )
                {
                    consVars[cell][k + 1] = density[0] * velocity[0][k];
                    for( int l = 0; l < 3; l++ )
                    {
                        consVars[cell][4 + 3 * k + l] = distortion[0](k, l);
                    }
                }
                consVars[cell][13] = density[0] * 
                    ( microEnergy( density[0], pressure[0] ) +
                      mesoEnergy( distortion[0] ) + 
                      macroEnergy( velocity[0] ) );
            }
            else
            {
                consVars[cell][0] = density[1];
                for( int k = 0; k < 3; k++ )
                {
                    consVars[cell][k + 1] = density[1] * velocity[1][k];
                    for( int l = 0; l < 3; l++ )
                    {
                        consVars[cell][4 + 3 * k + l] = distortion[1](k, l);
                    }
                }
                consVars[cell][13] = density[1] * 
                    ( microEnergy( density[1], pressure[1] ) +
                      mesoEnergy( distortion[1] ) + 
                      macroEnergy( velocity[1] ) );
            }
        }
    }
}

void HyperbolicPeshkovRomenski::boundaryConditions( BoundaryCondition type[4] )
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
            copyTo = cell; 

            switch (type[0] )
            {
                case transmissive: 
                    copyFrom = copyTo + M; 
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyFrom = copyTo + ( 2 * ( nGhostCells - i ) - 1 ) * M;
                    consVars[copyTo] = consVars[copyFrom];
                    consVars[copyTo][1] = - consVars[copyTo][1];
                    consVars[copyTo][2] = - consVars[copyTo][2];
                    break;
                case periodic: 
                    copyFrom = copyTo + M;
                    if( consVars[nCellsTot - M + j][1] > 0.0 )
                    {
                        if( i == 0 )
                            copyFrom = nCellsTot - M + j;
                        else
                            copyFrom = copyTo - M; 
                    }
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }

            // Right boundary (x = xMax)
            copyTo = nCellsTot - 1 - cell; 

            switch (type[1] )
            {
                case transmissive: 
                    copyFrom = copyTo - M; 
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyFrom = copyTo - ( 2 * ( nGhostCells - i ) - 1 ) * M;
                    consVars[copyTo] = consVars[copyFrom];
                    consVars[copyTo][1] = - consVars[copyTo][1];
                    consVars[copyTo][2] = - consVars[copyTo][2];
                    break;
                case periodic: 
                    copyFrom = copyTo - M;
                    if( consVars[M - 1 - j][1] < 0.0 )
                    {
                        if( i == 0 )
                            copyFrom = M - 1 - j;
                        else
                            copyFrom = copyTo + M; 
                    }
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
            copyTo = cell; 

            switch( type[2] )
            {
                case transmissive: 
                    copyFrom = copyTo + 1; 
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyFrom = copyTo + ( 2 * ( nGhostCells - j ) - 1 );
                    consVars[copyTo] = consVars[copyFrom];
                    consVars[copyTo][1] = - consVars[copyTo][1];
                    consVars[copyTo][2] = - consVars[copyTo][2];
                    break;
                case periodic: 
                    copyFrom = copyTo + 1;
                    if( consVars[( i + 1 ) * M - 1][2] > 0.0 )
                    {
                        if( j == 0 )
                            copyFrom = ( i + 1 ) * M - 1;
                        else
                            copyFrom = copyTo - 1; 
                    }
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }

            // Top boundary (y = yMax)
            copyTo = cell + M - 1 - 2 * j;

            switch( type[3] )
            {
                case transmissive: 
                    copyFrom = copyTo - 1; 
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case reflective:
                    copyFrom = copyTo - ( 2 * ( nGhostCells - j ) - 1 );
                    consVars[copyTo] = consVars[copyFrom];
                    consVars[copyTo][1] = 2.0 * 1.0 - consVars[copyTo][1];
                    consVars[copyTo][2] = - consVars[copyTo][2];
                    break;
                case periodic: 
                    copyFrom = copyTo - 1;
                    if( consVars[i * M][2] < 0.0 )
                    {
                        if( j == 0 )
                            copyFrom = copyTo + 1 - M; 
                        else
                            copyFrom = copyTo + 1;
                    }
                    consVars[copyTo] = consVars[copyFrom];
                    break;
                case constant:
                    break;
            }
        }
    }
}

void HyperbolicPeshkovRomenski::xSweep( double dt )
{
    std::vector< SimpleArray< double, 14 > > tempVars( nCellsTot );
#pragma omp parallel for 
    for( int i = 0; i < nCellsTot; i++ )
    {
        tempVars[i] = consVars[i];
    } 
    
    int cell; 
    SimpleArray< double, 14 > F_L, F_R;

#pragma omp parallel for private( cell, F_L, F_R )
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j; 
            slicFlux( dt, dx, 0, 
                    tempVars[cell - 2 * ( nCellsY + 2 * nGhostCells )],
                    tempVars[cell - ( nCellsY + 2 * nGhostCells )],
                    tempVars[cell], 
                    tempVars[cell + ( nCellsY + 2 * nGhostCells )], 
                    F_L );
            slicFlux( dt, dx, 0,
                    tempVars[cell - ( nCellsY + 2 * nGhostCells )],
                    tempVars[cell], 
                    tempVars[cell + ( nCellsY + 2 * nGhostCells )], 
                    tempVars[cell + 2 * ( nCellsY + 2 * nGhostCells )], 
                    F_R );
            consVars[cell] = tempVars[cell] + dt / dx * ( F_L - F_R ); 
        }
    } 
    renormalizeDistortion();
}

void HyperbolicPeshkovRomenski::ySweep( double dt )
{
    std::vector< SimpleArray< double, 14 > > tempVars( nCellsTot );
#pragma omp parallel for 
    for( int i = 0; i < nCellsTot; i++ )
    {
        tempVars[i] = consVars[i];
    } 
    
    int cell; 
    SimpleArray< double, 14 > F_B, F_T;
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
    renormalizeDistortion();
}

void HyperbolicPeshkovRomenski::addNonconservative( double dt )
{
    std::vector< SimpleArray< double, 14 > > tempVars( nCellsTot );
#pragma omp parallel for 
    for( int i = 0; i < nCellsTot; i++ )
    {
        tempVars[i] = consVars[i];
    } 
    
    int cell;
    SimpleArray< double, 14 > N;
#pragma omp parallel for private( cell, N )
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j; 
            // TODO: use boundary extrapolated values insted of FD for gradient
/*            finiteDiffNoncons( tempVars[cell], 
                    tempVars[cell - ( nCellsY + 2 * nGhostCells )], 
                    tempVars[cell + ( nCellsY + 2 * nGhostCells )], 
                    tempVars[cell - 1], 
                    tempVars[cell + 1], 
                    N ); */
            boundaryExtrapolatedNoncons( tempVars[cell], 
                    tempVars[cell - ( nCellsY + 2 * nGhostCells )], 
                    tempVars[cell + ( nCellsY + 2 * nGhostCells )], 
                    tempVars[cell - 1], 
                    tempVars[cell + 1], 
                    N ); 
            consVars[cell] = tempVars[cell] + dt * N; 
        }
    } 
    renormalizeDistortion();
}

void HyperbolicPeshkovRomenski::renormalizeDistortion()
{
    double rho;
    Eigen::Matrix3d A; 
    double scaleFactor; 
#pragma omp parallel for private( rho, A, scaleFactor )
    for( int i = 0; i < nCellsTot; i++ )
    {
        rho = getDensity( consVars[i] );
        A = getDistortion( consVars[i] );
        scaleFactor = pow( rho / ( rho_0 * A.determinant() ), 1.0 / 3.0 );
        A *= scaleFactor;
        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                consVars[i][4 + 3 * j + k] = A(j, k);
            }
        }
    }
}

void HyperbolicPeshkovRomenski::output2D( char* filename )
{
    int cell; 
    double x, y;
    double rho; 
    SimpleArray< double, 3 > u;
    double p; 
    Eigen::Matrix3d sigma;

    Eigen::Matrix3d curlyguy;

    std::ofstream fs; 
    fs.open( filename );
    fs << "x" << "\t" << "y" << "\t" << "rho" << "\t" << "p" << "\t" 
        << "u" << "\t" << "v" << "\t" << "w" << "\t" << "sigma11" << "\t" 
        << "sigma12" << "\t" << "sigma13" << "\t" << "sigma22" << "\t" 
        << "sigma23" << "\t" << "sigma33" << "\t" << "curlynorm" << std::endl; 

    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j;
            x = domain[0] + (i - nGhostCells + 0.5 ) * dx; 
            y = domain[2] + (j - nGhostCells + 0.5 ) * dy; 
            rho = getDensity( consVars[cell] );
            u = getVelocity( consVars[cell] );
            p = getPressure( consVars[cell] );
            sigma = getShearStress( consVars[cell] );
            curlyguy = getCurlTerm( consVars[cell], 
                consVars[cell - nCellsY - 2 * nGhostCells], 
                consVars[cell + nCellsY + 2 * nGhostCells], 
                consVars[cell - 1], 
                consVars[cell + 1] );

            fs << x << " \t" << y << " \t" 
                << rho << " \t" << p << " \t" << u[0] << " \t" << u[1] << " \t" 
                << u[2] << " \t" << sigma(0, 0) - p << "\t" << sigma(0, 1) << "\t"
                << sigma(0, 2) << "\t" << sigma(1, 1) - p << "\t" 
                << sigma(1, 2) << "\t" << sigma(2, 2) - p << "\t" 
                << curlyguy.norm() << std::endl;
        }
        fs << std::endl; 
    }
    fs.close();
}

void HyperbolicPeshkovRomenski::output1DSliceX( char* filename )
{
    int cell;
    double x;
    double rho; 
    SimpleArray< double, 3 > u;
    double uAbs;
    double p; 
    double e; 

    std::ofstream fs; 
    fs.open( filename );
    fs << "x" << "\t" << "rho" << "\t" << "u" << "\t" << "v" << "\t" 
        << "p" << "\t" << "e" << "\t" << "uMax" << std::endl;

    int j = nGhostCells + nCellsY / 2;
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        cell = i * ( nCellsY + 2 * nGhostCells ) + j;
        x = domain[0] + (i - nGhostCells + 0.5 ) * dx; 
        rho = getDensity( consVars[cell] );
        u = getVelocity( consVars[cell] );
        uAbs = sqrt( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
        p = getPressure( consVars[cell] );
        e = microEnergy( rho, p );

        fs << x << "\t" << rho << "\t" << u[0] << "\t" << u[1] << "\t" 
            << p << "\t" << e << "\t" << uAbs << std::endl;
    }
    fs.close();
}

void HyperbolicPeshkovRomenski::output1DSliceY( char* filename )
{
    int cell;
    double y;
    double rho; 
    SimpleArray< double, 3 > u;
    double uAbs;
    double p; 
    double e; 

    std::ofstream fs; 
    fs.open( filename );
    fs << "y" << "\t" << "rho" << "\t" << "u" << "\t" << "v" << "\t" 
        << "p" << "\t" << "e" << "\t" << "uAbs" << std::endl;

    int i = nGhostCells + nCellsX / 2;
    for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
    {
        cell = i * ( nCellsY + 2 * nGhostCells ) + j;
        y = domain[2] + (j - nGhostCells + 0.5 ) * dy; 
        rho = getDensity( consVars[cell] );
        u = getVelocity( consVars[cell] );
        uAbs = sqrt( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
        p = getPressure( consVars[cell] );
        e = microEnergy( rho, p );

        fs << y << "\t" << rho << "\t" << u[0] << "\t" << u[1] << "\t" 
            << p << "\t" << e << "\t" << uAbs << std::endl;
    }
    fs.close();
}

/* CLASS HPR_FLUID */

HPR_Fluid::HPR_Fluid( double _shearSoundSpeed, double _referenceDensity, 
        int _nCellsX, int _nCellsY, double _domain[4], 
        double _gamma, double _strainDissipationTime ) : 
    HyperbolicPeshkovRomenski( _shearSoundSpeed, _referenceDensity, 
            _nCellsX, _nCellsY, _domain) 
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

double HPR_Fluid::microEnergy( double density, double pressure )
{
    return pressure / ( ( gamma - 1.0 ) * density );
}

double HPR_Fluid::getTimeStep( const double c_CFL )
{ 
    std::vector< double > Sx( nCellsTot ); 
    std::vector< double > Sy( nCellsTot ); 
    double rho;
    SimpleArray< double, 3 > u; 
    double p;
    double a; 

#pragma omp parallel for private( rho, u, p, a )
    // TODO: reduce for loop with max
    for( int i = 0; i < nCellsTot; i++ )
    {
        rho = getDensity( consVars[i] );
        u = getVelocity( consVars[i] );
        p = getPressure( consVars[i] );
        a = sqrt( gamma * p / rho + 4.0 * c_s * c_s / 3.0 );
        Sx[i] = fabs( u[0] ) + a; 
        Sy[i] = fabs( u[1] ) + a; 
    }
    double SxMax = *std::max_element( Sx.begin(), Sx.end() );
    double SyMax = *std::max_element( Sy.begin(), Sy.end() );
    return c_CFL * std::min( dx / SxMax, dy / SyMax);
}

void HPR_Fluid::integrateODE( double dt )
{
    int cell; 
    double tol = 1.0e-16;

#pragma omp parallel for private( cell )
    for( int i = nGhostCells; i < nGhostCells + nCellsX; i++ )
    {
        for( int j = nGhostCells; j < nGhostCells + nCellsY; j++ )
        {
            cell = i * ( nCellsY + 2 * nGhostCells ) + j; 
            integrate_adaptive( make_controlled( tol, tol, stepper_type() ),
                    HPR_ODE( tau, rho_0 ), consVars[cell], 0.0, 0.0 + dt, 
                    1.0e-6 * dt );
        }
    }
    renormalizeDistortion();
}

/* CLASS HPR_SOLID */

HPR_Solid::HPR_Solid( double _shearSoundSpeed, double _referenceDensity, 
        int _nCellsX, int _nCellsY, double _domain[4], 
        double _c_0, double _Gamma_0, double _s_H ) : 
    HyperbolicPeshkovRomenski( _shearSoundSpeed, _referenceDensity, 
            _nCellsX, _nCellsY, _domain) 
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
{ // TODO: FIX FOR 2D!
    return 0.0;
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
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            A(i, j) = Q[4 + 3 * i + j];
        }
    }

    Eigen::Matrix3d G = A.transpose() * A;
    Eigen::Matrix3d devG = G - G.trace() / 3.0 * Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Psi = - 3.0 * pow( nu, 5.0 / 3.0 ) / tau * A * devG;

    // TODO: std::fill to set elements of S to zero by default
    for( int i = 0; i < 4; i++ )
    { // no source terms for density, velocity
        S[i] = 0.0;
    }

    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
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

void configurate( const char* inputFile, int& nCellsX, int& nCellsY,
        double& CFL, double& tStop, double& c_s, double& rho_0, double& gamma,
        double& tau, double domain[4], double& initDiscontPos, 
        Direction& initDiscontDir, BoundaryCondition BCs[4], 
        double initDensity[2], SimpleArray< double, 3 > initVelocity[2],
        Eigen::Matrix3d initDistortion[2], double initPressure[2] )
{
    // Read input file
    libconfig::Config cfg; 
    cfg.readFile( inputFile );

    double mu;
    int BCs_int[4], dir_int; 

    cfg.lookupValue( "nCellsX", nCellsX );
    cfg.lookupValue( "nCellsY", nCellsY );
    cfg.lookupValue( "CFL", CFL );
    cfg.lookupValue( "tStop", tStop );
    cfg.lookupValue( "shearWaveSpeed", c_s );
    cfg.lookupValue( "referenceDensity", rho_0 );
    cfg.lookupValue( "heatCapacityRatio", gamma );
    cfg.lookupValue( "viscosity", mu );
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
    cfg.lookupValue( "w_1", initVelocity[0][2] );
    cfg.lookupValue( "p_1", initPressure[0] );
    cfg.lookupValue( "rho_2", initDensity[1] );
    cfg.lookupValue( "u_2", initVelocity[1][0] );
    cfg.lookupValue( "v_2", initVelocity[1][1] );
    cfg.lookupValue( "w_2", initVelocity[1][2] );
    cfg.lookupValue( "p_2", initPressure[1] );

    for( int i = 0; i < 4; i++ )
    {
        BCs[i] = (BoundaryCondition) BCs_int[i];
    }

    initDiscontDir = (Direction) dir_int;

    tau = 6.0 * mu / (rho_0 * c_s * c_s );

    for( int i = 0; i < 2; i++ )
    {
        initDistortion[i] = pow( initDensity[i] / rho_0, 1.0 / 3.0 ) *
            Eigen::Matrix3d::Identity();
    }
}

#endif
