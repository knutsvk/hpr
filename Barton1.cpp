#include "HPR.h"

using namespace std;

int main( int argc, char* argv[] )
{
    int N = 200;
    if( argc > 1 ) N = atoi( argv[1] );
    double dom[2] = { 0.00, 0.01 };
    double c = 0.9;

    double x0 = ( dom[1] - dom[0] ) / 2.0;

    double rho_L = 8930.0;
    double rho_R = 8930.0;

    double u_L[3] = { 0.0, 0.0, 0.0 };
    double u_R[3] = { 0.0, 0.0, 0.0 };

    Matrix3d A_L, A_R; 
    A_L << 0.95, 0.00, 0.00,
        0.00, 0.00, 0.00, 
        0.00, 0.00, 1.00;
    A_R << 1.00, 0.00, 0.00,
        0.00, 1.00, 0.00, 
        0.00, 0.00, 1.00;

    double p_L = ;
    double p_R = ;

    double tStop = 0.06e-6;

    // Copper
    double C_s = 2190.0; // Barton: 2100
    double rho0 = 8930.0; // Barton: 8930, Wikipedia: 8960
    double C_0 = 3940.0 ; // Barton: 4600, Wikipedia: 3933
    double G0 = 2.00; // Wikipedia: 1.99
    double s = 1.48; // Wikipedia: 1.50

    HPR_Solid state( C_s, rho0, N, dom, C_0, G0, s );
    state.initialize( x0, rho_L, rho_R, u_L, u_R, A_L, A_R, p_L, p_R );

    double t = 0.0;
    double dt;
    int iter = 1;
    while(t < tStop)
    {
        dt = state.timeStep(c);

        if(iter < 10)
        {
            dt *= 0.1 / c;
        }

        if(t + dt > tStop)
        {
            dt = tStop - t;
        }

        state.transmissiveBCs();
        state.slic(dt);
        state.advancePDE(dt);

        t += dt;
        iter++;
    }

    state.output();
}
