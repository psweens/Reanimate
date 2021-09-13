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


// Approx. cortical pressures for diameters < ~90um (based on Lorthois et al. 2011)
void Vasculature::cortexBoundaryPress() {

    int seg{};
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 0) {
            seg = nodseg(0,bcnod(inodbc));
            if (BCgeo(inodbc) == 1) {
                bcprfl(inodbc) = -0.0287*pow(diam(seg),2) + 2.84884*diam(seg) + 19.966;
                bctyp(inodbc) = 0;
            }
            else if (BCgeo(inodbc) == 3) {
                bcprfl(inodbc) = -1.e-4 * pow(diam(seg), 3) + 0.0186 * pow(diam(seg),2) - 1.0616 * diam(seg) + 35.04;
                bctyp(inodbc) = 0;
            }
            cout<<bcnodname(inodbc)<<"\t"<<bcprfl(inodbc)<<"\t"<<diam(seg)<<endl;
        }
    }

}