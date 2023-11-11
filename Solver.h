#ifndef MINREP_SOLVER_H
#define MINREP_SOLVER_H

#include <Eigen/Eigen>
#include "Body.h"

using namespace std;
using namespace Eigen;

class Solver {

private:

    const double G = 6.674 * pow(10, -11);
    const int subSteps = 16;

    vector<Body> bodyList;

public:

    Solver() {
        this->bodyList = initialConditions();
    }

    void passTime(double dt) {

        for (int i = 0; i < subSteps; i++) {
            doVerlet(dt/subSteps);
        }

    }

    //executes verlet integration on all bodies
    void doVerlet(double dt) {

        for (Body &body : bodyList)
            body.acceleration = totalGravitationalAccelerationOf(body);


        for (Body &body : bodyList)
            body.verletNextStep(dt);

    }

    //calculates potential energy
    long double calcPot() {

        long double total = 0;

        for (Body &body1 : bodyList) {
            for (Body &body2: bodyList) {
                if (&body1 != &body2) {
                    total += ((long double)(-G * body1.mass * body2.mass)) / (long double)body1.vectorTo(body2).norm();
                }
            }
        }

        return total;

    }

    Vector2<long double> totalGravitationalAccelerationOf(Body &body) {

        Vector2<long double> acc(0, 0);

        for (Body &otherBody : bodyList)
            if (body.vectorTo(otherBody).norm() != 0)
                acc += gravitationalAccelerationOf(body, otherBody);

        return acc;
    };

    //acceleration exerted by body2, on body1
    Vector2<long double> gravitationalAccelerationOf(Body &body1, Body &body2) {

        Vector2<long double> rHat = body1.vectorTo(body2).normalized();
        long double rSquared = body1.vectorTo(body2).squaredNorm();

        return rHat * G * body2.mass/rSquared;
    };

    //calculates important quantities
    vector<long double> quantsInfo() {

        long double momx = 0;
        long double momy = 0;
        long double kin = 0;
        long double energy = 0;

        long double pot = calcPot();

        for (Body body : bodyList) {
            momx += body.momentum().x();
            momy += body.momentum().y();
            kin += body.kineticEnergy();
        }

        energy = kin + pot;

        return {momx, momy, energy, kin, pot};
    };

    //returns initial conditions of system
    vector<Body> initialConditions() {

        vector<Body> list;

        Vector2<long double> pos = Vector2<long double>(-1000, 1000);
        Vector2<long double> vel = Vector2<long double>(50, -100);


        list.emplace_back(70000000000000000.0,
                          pos,
                          vel);

        pos = Vector2<long double>(1000, 1000);
        vel = Vector2<long double>(-100, -50);


        list.emplace_back(50000000000000000.0,
                          pos,
                          vel);

        pos = Vector2<long double>(1000, -1000);
        vel = Vector2<long double>(-50, 100);

        list.emplace_back(50000000000000000.0,
                          pos,
                          vel);

        return list;

    }

};


#endif
