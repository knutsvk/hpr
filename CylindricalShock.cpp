#include "HPR.h"
#include <libconfig.h++>

using namespace Eigen; 
using namespace std;

int main( int argc, char* argv[] )
{
    libconfig::Config cfg; 
    cfg.readFile("CylindricalShock.cfg");
    // Computational parameters
    int N;
    double c;
    double dom[4];
    cfg.lookupValue( "nCellsPerDirection", N );
    cfg.lookupValue( "CFL", c );
    cfg.lookupValue( "xMin", dom[0] );
    cfg.lookupValue( "xMax", dom[1] );
    cfg.lookupValue( "yMin", dom[2] );
    cfg.lookupValue( "yMax", dom[3] );

    // Cylindrical shock tube 
    BoundaryCondition BCs[4] = { transmissive, transmissive, transmissive,
        transmissive };
    Direction dir;
    cfg.lookupValue( "initialDiscontinuityDirection", dir );
    double tStop = 0.25;
    double R = 0.4;
    double rho_L = 1.000;
    double rho_R = 0.125;
    SimpleArray< double, 3 > u_L = { 0.0, 0.0, 0.0 };
    SimpleArray< double, 3 > u_R = { 0.0, 0.0, 0.0 };
    Matrix3d A_L = Matrix3d::Identity();
    Matrix3d A_R = 0.5 * Matrix3d::Identity();
    double p_L = 1.0;
    double p_R = 0.1;

    // Air properties
    double C_s = 0.25; 
    double rho0 = 1.0;
    double g = 1.4;
    double t_PSL = 6.0e-4 / ( rho0 * C_s * C_s );

    HPR_Fluid state( C_s, rho0, N, N, dom, g, t_PSL );
    state.initialize( R, dir, rho_L, rho_R, u_L, u_R, A_L, A_R, p_L, p_R );

    double t = 0.0;
    double dt;
    int iter = 1;
    while(t < tStop)
    {
        dt = state.getTimeStep( c );
        if( iter < 10 )
            dt *= 0.1 / c;
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

    char filename[50]; 
    sprintf( filename, "CylindricalSlices_N%d.out", N );
    state.output1DSlices( filename );
    sprintf( filename, "CylindricalSurface_N%d.out", N );
    state.output2D( filename );
}
