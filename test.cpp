#include "Materials.h"

using namespace std;

int main(int argc, char* argv[])
{
    int N = 5000;
    if(argc > 1) N = atoi(argv[1]);
    double dom[2] = {0.0, 1.0};
    double c = 0.9;

    double x0 = 0.5;
    double rho[2] = {1.0, 0.125};
    double u[2] = {0.0, 0.0};
    double p[2] = {1.0, 0.1};

    double tStop = 0.25;

    Material state(N, dom);
    state.initialize(x0, rho, u, p);

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
        state.force(dt);
        state.advancePDE(dt);

        t += dt;
        iter++;
    }

    state.output();
}
