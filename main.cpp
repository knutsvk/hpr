#include "HPR.h"

using namespace Eigen; 
using namespace std;

int main( int argc, char* argv[] )
{
    // Declaration of variables
    int Nx, Ny, iter;
    double c, tStop, C_s, rho0, g, t_PSL, dom[4], x_0, rho[2], p[2], t, dt, l;
    SimpleArray< double, 3 > u[2];
    Matrix3d A[2];
    BoundaryCondition BCs[4];
    Direction dir; 
    char sim[100], infile[100], outfile[100];

    // Find which simulation to run
    if( argc == 2 )
        sprintf( sim, "%s", argv[1] );
    else
    {
        cout << "Available simulations: " << endl << endl
            << "CylindricalShock" << endl 
            << "DoubleShearLayer" << endl
            << "LaminarBoundaryLayer" << endl
            << "LidDrivenCavity" << endl
            << "PeriodicWave" << endl
            << "StokesFirstProblem" << endl
            << "Toro1" << endl
            << "Toro2" << endl
            << "Toro3" << endl
            << "Toro4" << endl 
            << endl << "Enter name of simulation to run: ";
        cin >> sim; 
    }
    sprintf( infile, "%s.cfg", sim );

    // Get values from configuration file
    configurate( infile, Nx, Ny, c, tStop, C_s, rho0, g, t_PSL, dom, x_0,
            dir, BCs, rho, u, A, p );

    // Allocate memory, set material properties
    HPR_Fluid state( C_s, rho0, Nx, Ny, dom, g, t_PSL );

    // Apply initial conditions
    if( !strcmp( sim, "DoubleShearLayer" ) )
        state.initialiseDoubleShearLayer();
    else if( !strcmp( sim, "ConvergenceStudies" ) )
    {
        state.exactConvergenceSolution();
        state.initialiseConvergenceTest();
    }
    else
        state.initialise( x_0, dir, rho, u, A, p );

    // Enter simulation loop
    t = 0.0;
    iter = 1;
    l = 0.0;
    while(t < tStop)
    {
        cout << "iter: " << iter << ", t = " << t;
        
        if( t / tStop >= l / 100.0 )
        {
            sprintf( outfile, "./Results/%s_1DX_Nx%d_Ny%d_%d.out", sim, Nx, Ny, (int) l );
            state.output1DSliceX( outfile );
            sprintf( outfile, "./Results/%s_1DY_Nx%d_Ny%d_%d.out", sim, Nx, Ny, (int) l );
            state.output1DSliceY( outfile );
            sprintf( outfile, "./Results/%s_2D_Nx%d_Ny%d_%d.out", sim, Nx, Ny, (int) l );
            state.output2D( outfile );
            l += 10.0;
        }

        dt = state.getTimeStep( c );
        if( iter < 10 )
            dt *= 0.1;
        if( t + dt > tStop )
            dt = tStop - t;

        cout << ", dt = " << dt;

/*        state.boundaryConditions( BCs );
        state.integrateODE( 0.5 * dt );
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }

        state.boundaryConditions( BCs );
        state.xSweep( 0.5 * dt );
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }

        state.boundaryConditions( BCs );
        state.ySweep( dt );
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }

        state.boundaryConditions( BCs );
        state.xSweep( 0.5 * dt );
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }

        state.boundaryConditions( BCs );
        state.integrateODE( 0.5 * dt );
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }

        state.boundaryConditions( BCs );
        state.diffuse();
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }

        state.renormalizeDistortion();
        if( !state.isPhysical() )
        {
            cout << "Unphysical state encountered in iteration " << iter 
                << ", time = " << t << ". " << endl;
            return 1;
        }*/

        state.boundaryConditions( BCs );
        state.xSweep( dt );
        state.ySweep( dt );
        state.integrateODE( dt );
        state.renormalizeDistortion();
        state.diffuse();

        cout << endl; 

        t += dt;
        iter++;
    }
    cout << "Done!" << endl; 

    // Write final results to file
    sprintf( outfile, "./Results/%s_1DX_Nx%d_Ny%d_%d.out", sim, Nx, Ny, (int) l );
    state.output1DSliceX( outfile );
    sprintf( outfile, "./Results/%s_1DY_Nx%d_Ny%d_%d.out", sim, Nx, Ny, (int) l );
    state.output1DSliceY( outfile );
    sprintf( outfile, "./Results/%s_2D_Nx%d_Ny%d_%d.out", sim, Nx, Ny, (int) l );
    state.output2D( outfile );

    return 0;
}
