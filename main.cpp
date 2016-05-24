#include "HPR.h"

using namespace Eigen; 
using namespace std;

int main( int argc, char* argv[] )
{
    // Declaration of variables
    int Nx, Ny, iter;
    double c, tStop, C_s, rho0, g, t_PSL, dom[4], x_0, rho[2], p[2], t, dt;
    SimpleArray< double, 3 > u[2];
    Matrix3d A[2];
    BoundaryCondition BCs[4];
    Direction dir; 
    char filename[50], inputFile[50], outputFile[50];

    // Find which simulation to run
    if( argc == 2 )
        sprintf( filename, "%s", argv[1] );
    else
    {
        cout << "Available simulations: " << std::endl 
            << "CylindricalShock" << std::endl 
            << "PeriodicWave" << std::endl
            << "StokesFirstProblem" << std::endl
            << "Enter name of simulation to run: ";
        cin >> filename; 
    }
    sprintf( inputFile, "%s.cfg", filename );

    // Get values from configuration file
    configurate( inputFile, Nx, Ny, c, tStop, C_s, rho0, g, t_PSL, dom, x_0,
            dir, BCs, rho, u, A, p );

    // Allocate memory, set material properties
    HPR_Fluid state( C_s, rho0, Nx, Ny, dom, g, t_PSL );

    // Apply initial conditions
    state.initialize( x_0, dir, rho, u, A, p );

    // Enter simulation loop
    t = 0.0;
    iter = 1;
    while(t < tStop)
    {
        dt = state.getTimeStep( c );
        if( iter < 10 )
            dt *= 0.1;
        if( t + dt > tStop )
            dt = tStop - t;

        state.boundaryConditions( BCs );

        state.integrateODE( 0.5 * dt );
        state.xSweep( 0.5 * dt );
        state.ySweep( dt );
        state.xSweep( 0.5 * dt );
        state.integrateODE( 0.5 * dt );

        t += dt;
        iter++;
    }

    // Write results to file
    sprintf( outputFile, "%s_1DSliceX_Nx%d_Ny%d.out", filename, Nx, Ny );
    state.output1DSliceX( outputFile );
    sprintf( outputFile, "%s_1DSliceY_Nx%d_Ny%d.out", filename, Nx, Ny );
    state.output1DSliceY( outputFile );
    sprintf( outputFile, "%s_2D_Nx%d_Ny%d.out", filename, Nx, Ny );
    state.output2D( outputFile );
}
