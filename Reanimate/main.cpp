#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"

#include "omp.h"

#include "MicroCell.hpp"
#include "DiscreteContinuum.hpp"
#include <sys/resource.h>


using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    // Generate discrete flow solution
    DiscreteContinuum hybrid;
    hybrid.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate-DC/Build_Data/";
    hybrid.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate-DC/Load_Data/";
    hybrid.discreteNet.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate-DC/Build_Data/";
    hybrid.discreteNet.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate-DC/Load_Data/";
    hybrid.discreteNet.setBuildPath(true);
    hybrid.discreteNet.loadNetwork("1Network.dat");
    hybrid.discreteNet.bloodFlow(true);
    hybrid.discreteNet.pictureNetwork("Network_Pressure.ps",hybrid.discreteNet.segpress);

    hybrid.discreteNet.findBoundingBox();

    // Generate micro-cell
    hybrid.cell.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate-DC/Build_Data/";
    hybrid.cell.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate-DC/Load_Data/";
    hybrid.cell.setBuildPath(false); // Create build directory / delete contents

    hybrid.cell.setDiamDistrib(hybrid.discreteNet.diam(find(hybrid.discreteNet.vesstyp == 2)));
    hybrid.cell.setLengthDistrib(hybrid.discreteNet.lseg(find(hybrid.discreteNet.vesstyp == 2)));
    hybrid.cell.rotationAngle = 0.;
    hybrid.cell.computeConductivity("crossCell3D");

    // Run vessel classification
    imat inOutlets = zeros<imat>(2,2);
    inOutlets(0,0) = 830;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 825;
    inOutlets(1,1) = 2;
    //inOutlets(2,0) = 824;
    //inOutlets(2,1) = 1;
/*    inOutlets(0,0) = 78;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 76;
    inOutlets(1,1) = 2;
    inOutlets(2,0) = 225;
    inOutlets(2,1) = 1;*/

    // Generate spatial graph
    hybrid.graph.generate(hybrid.discreteNet, true);
    hybrid.graph.analyseTopology(inOutlets);

    // Run discrete-continuum model
    hybrid.runHybrid();

}
