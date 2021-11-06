#include "Vasculature.hpp"

using namespace reanimate;

void Vasculature::assignBoundaryHD() {

    double consthd = 0.4;
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(0., 1.0);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++)    {
        if (bctyp(inodbc) == 3) {
            double p = dis(gen);
            if (p >= 0 && p <= consthd) {bchd(inodbc) = p / (3 * pow(consthd,2));}
            else if (p > consthd && p < 1.5*consthd)    {bchd(inodbc) = 2 * (3*consthd / 2-p) / (3 * pow(consthd,2));}
            else {bchd(inodbc) = 0.;}
        }
    }

}
