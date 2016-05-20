#include "HPR.h"

using namespace Eigen; 
using namespace std;

int main( int argc, char* argv[] )
{
    bool printProgress = false; 

    int N = 200;
    if( argc > 1 ) N = atoi( argv[1] );
    double dom[2] = { 0.0, 1.0 };
    double c = 0.9;

    double x0 = ( dom[1] - dom[0] ) / 2.0;

    double rho_L = 1.000;
    double rho_R = 0.125;
    SimpleArray< double, 3 > u_L = { 0.0, 0.0, 0.0 };
    SimpleArray< double, 3 > u_R = { 0.0, 0.0, 0.0 };

    Matrix3d A_L = Matrix3d::Identity();
    Matrix3d A_R = 0.5 * Matrix3d::Identity();

    double p_L = 1.0;
    double p_R = 0.1;

    double tStop = 0.25;

    // Air
    double C_s = 0.5; 
    double rho0 = 1.0;
    double g = 1.4;
    double t_PSL = 6.0e-4 / ( rho0 * C_s * C_s );

    if( printProgress )
        cout << "Toro's test 1 in the HPR framework with " << N << " cells. " << endl;

    if( printProgress )
        cout << "Allocating memory... " << endl; 

    HPR_Fluid state( C_s, rho0, N, dom, g, t_PSL );

    if( printProgress )
        cout << "Initializing test case... " << endl; 

    state.initialize( x0, rho_L, rho_R, u_L, u_R, A_L, A_R, p_L, p_R );

    if( printProgress )
        cout << "Simulating... " << endl << "iter \t time" << endl;

    double t = 0.0;
    double dt;
    int iter = 1;
    while(t < tStop)
    {
        if( printProgress )
            cout << iter << " \t " << t << endl; 

        dt = state.getTimeStep( c );

        if( iter < 10 )
            dt *= 0.1 / c;

        if( t + dt > tStop )
            dt = tStop - t;

        state.transmissiveBCs();

        state.integrateODE( 0.5 * dt );

        state.slic( dt );

        state.advancePDE( dt );

        state.integrateODE( 0.5 * dt );

        t += dt;
        iter++;
    }

    state.output();

    if( printProgress )
        cout << "Done!" << endl; 
}
