#include "HPR.h"
#include <cstring>

using namespace std;

int main( int argc, char* argv[] )
{
    // Declaration of variables
    int Nx, Ny, iter;
    double c, tStop, g, dom[4], x_0, rho[2], p[2], t, dt;
    SimpleArray< double, 2 > u[2];
    BoundaryCondition BCs[4];
    Direction dir;
    char sim[100], fin[100], fout[100];

    // Find which simulation to run
    if( argc == 2 )
        sprintf( sim, "%s", argv[1] );
    else
    {
        cout << "Available simulations: " << endl << endl
            << "ConvergenceStudies" << endl
            << "CylindricalShock" << endl
            << "PeriodicWave" << endl
            << "Toro1" << endl
            << "Toro2" << endl
            << "Toro3" << endl
            << "Toro4" << endl
            << endl << "Enter name of simulation to run: ";
        cin >> sim;
    }
    sprintf( fin, "%s.cfg", sim );

    // Get values from configuration file
    configurate( fin, Nx, Ny, c, tStop, g, dom, x_0, dir, BCs, rho, u, p );

    // Allocate memory, set material properties
    Euler state( Nx, Ny, dom, g );

    // Apply initial conditions
    if( !strcmp( sim, "ConvergenceStudies" ) )
    {
        state.exactConvergenceSolution();
        state.initialiseConvergenceTest();
    }
    else
        state.initialise( x_0, dir, rho, u, p );

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
        state.xSweep( 0.5 * dt );
        state.ySweep( dt );
        state.xSweep( 0.5 * dt );

        cout << "iter: " << iter << ", t = " << t << ", dt = " << dt << endl;

        t += dt;
        iter++;
    }
    cout << "Done!" << endl;

    // Write final results to file
    sprintf( fout, "./Results/%s_1DX_Nx%d_Ny%d.out", sim, Nx, Ny );
    state.output1DSliceX( fout );
    sprintf( fout, "./Results/%s_1DY_Nx%d_Ny%d.out", sim, Nx, Ny );
    state.output1DSliceY( fout );
    sprintf( fout, "./Results/%s_2D_Nx%d_Ny%d.out", sim, Nx, Ny );
    state.output2D( fout );

    return 0;
}
