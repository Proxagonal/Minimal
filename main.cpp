#include "Solver.h"

using namespace std;
using namespace Eigen;


long double deviation(long double x, long double y);
void compare(vector<long double> now, vector<long double> init);

vector<string> names = {"Momentum_x: ", "Momentum_y: ", "Total Energy: ", "Kinetic: ", "Potential: "};

int main() {

    Solver solver;

    double dt = 0.001;
    int T = 10;

    vector<long double> initQuants = solver.quantsInfo();

    for (int i = 0; i < T/dt; i++) {

        solver.passTime(dt);

        //if an integer time has passed
        if (fmod(i*dt, 1) == 0) {

            vector<long double> quants = solver.quantsInfo();

            cout << "----------" << endl;
            cout << "TIME: " << i*dt << endl;
            compare(quants, initQuants);
        }

    }

    return 0;
}

//compares conserved quantities initially to now
void compare(vector<long double> now, vector<long double> init) {
    // kinetic energy and potential energy aren't supposed to be conserved, so I don't print them
    for (int i = 0; i<3; i++)
        cout << names.at(i) << deviation(now.at(i), init.at(i)) << endl;
}

long double deviation(long double x, long double y) {
    return (x-y)/y;
}

