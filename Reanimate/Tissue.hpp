//
// Created by sweene01 on 07/07/2020.
//

#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

#ifndef Tissue_hpp
#define Tissue_hpp

namespace reanimate {

    class Tissue {

    private:

        int nt;
        double volDt,accubgtt;
        vec vseg,r0,sumbgtv;
        mat snode,gfac;

        double bicgstab(vec &b, vec &x, int &n, double &eps, int &itmax);

    };

}

#endif /* Tissue_hpp */
