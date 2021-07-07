#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"

using namespace reanimate;
using namespace std;

int main(int argc, char** argv) {

    // Generate discrete flow solution
    DiscreteContinuum hybrid;
    hybrid.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.discreteNet.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.discreteNet.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.discreteNet.setBuildPath(true);
    hybrid.discreteNet.loadNetwork("Medulla.txt");

    hybrid.discreteNet.pictureNetwork("Network_Diameter.ps",hybrid.discreteNet.diam);
    mat extraD = zeros<mat>(hybrid.discreteNet.getNseg(), 1);
    hybrid.discreteNet.printAmira("amiraNetwork.am", extraD);

    hybrid.discreteNet.bloodFlow(true, false);
    hybrid.discreteNet.printVisuals();

    hybrid.discreteNet.findBoundingBox();

    // Run vessel classification
    // Mes. 1
    imat inOutlets = zeros<imat>(2,2);
    inOutlets(0,0) = 830;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 825;
    inOutlets(1,1) = 2;

    // Mes. 3
    /*inOutlets = zeros<imat>(2,2);
    inOutlets(0,0) = 78;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 76;
    inOutlets(1,1) = 2;*/

    // Mes. 5
    /*imat inOutlets = zeros<imat>(14,2);
    inOutlets(0,0) = 999;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 2996;
    inOutlets(1,1) = 2;
    inOutlets(2,0) = 3;
    inOutlets(2,1) = 1;
    inOutlets(3,0) = 994;
    inOutlets(3,1) = 2;
    inOutlets(4,0) = 19;
    inOutlets(4,1) = 2;
    inOutlets(5,0) = 992;
    inOutlets(5,1) = 2;
    inOutlets(6,0) = 1109;
    inOutlets(6,1) = 2;
    inOutlets(7,0) = 1102;
    inOutlets(7,1) = 2;
    inOutlets(8,0) = 579;
    inOutlets(8,1) = 1;
    inOutlets(9,0) = 984;
    inOutlets(9,1) = 2;
    inOutlets(10,0) = 992;
    inOutlets(10,1) = 1;
    inOutlets(11,0) = 550;
    inOutlets(11,1) = 2;
    inOutlets(12,0) = 456;
    inOutlets(12,1) = 1;
    inOutlets(13,0) = 297;
    inOutlets(13,1) = 2;*/

    // Cortex
    /*imat inOutlets = zeros<imat>(4,2);
    inOutlets(0,0) = 601;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 599;
    inOutlets(1,1) = 1;
    inOutlets(2,0) = 595;
    inOutlets(2,1) = 2;
    inOutlets(3,0) = 576;
    inOutlets(3,1) = 2;*/


    // Generate spatial graph
    hybrid.graph.generate(hybrid.discreteNet, true);
    //hybrid.graph.findTree(inOutlets);
    //hybrid.graph.mapClassification(hybrid.discreteNet);
    hybrid.graph.analyseTopology(inOutlets, hybrid.discreteNet);
    hybrid.discreteNet.analyseBoundaryType();
    hybrid.discreteNet.pictureNetwork("Network_Geometry.ps",conv_to<vec>::from(hybrid.discreteNet.vesstyp));


    // Generate micro-cell
    hybrid.cell.buildPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.cell.loadPath = "/home/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.cell.setBuildPath(false); // Create build directory / delete contents

    hybrid.cell.setDiamDistrib(hybrid.graph.diam(find(hybrid.graph.vesstyp == 2)));
    hybrid.cell.setLengthDistrib(hybrid.graph.lseg(find(hybrid.graph.vesstyp == 2)));
    hybrid.cell.rotationAngle = 0.;
    hybrid.cell.computeConductivity("crossCell3D");


    // Run discrete-continuum model
    hybrid.runHybrid();



}
