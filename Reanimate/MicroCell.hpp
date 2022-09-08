//
// Created by Paul Sweeney on 14/02/2021.
//

#include <string>
#include <armadillo>
#include <sys/resource.h>
#include "Vasculature.hpp"

using namespace arma;
using namespace std;

#ifndef MicroCell_hpp
#define MicroCell_hpp

namespace reanimate {

    class MicroCell : public Vasculature {

    public:

        bool cell2D;
        double omega{},kappa{},rotationAngle{},netVol{},vascDens{},lsegDens{},surfDens{},surfVolRatio{},R{},aniScaleY{1.},aniScaleZ{1.};
        vec diamDistrib, lengthDistrib, eucLengths, gridSpacing;
        mat conductivity,cellSegpress;
        void hexCell2D();
        void crossCell2D();
        void crossCell3D();

        void flowMicroCell();
        void computeConductivity(const string cellType, const int iterations=1e3);

        // Setter functions
        void setEdgeDiamDistrib(const vec d);
        void setEdgeLengthDistrib(const vec l);
        void setEucLengthDistrib(const vec l);
        void setRotAngle(const double angle);

        MicroCell();
        ~MicroCell();

    private:

        ivec Bin, Bout, bcPairs;
        vec vessOrient,tissForce;
        mat cA,cB,cC,cE,cF,BA,EA,unitCL;

        void setup_microCellArrays(bool tessellate=false);
        void setup_mcFlowArrays();
        void loadMicroCell();
        void findBCpairs();
        void assignGeometry();
        void analyseMicroCell();
        void printCellAnalysis(string filename);

    };

}

#endif /* MicroCell.hpp */
