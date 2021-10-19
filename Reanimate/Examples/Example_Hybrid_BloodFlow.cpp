#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"
#include "math.h"

using namespace reanimate;
using namespace std;

void example_Hybrid_BloodFlow() {

    // Set number of threads for any multithreading
    omp_set_num_threads(48);

    // Define discrete-continuum (hybrid) model
    DiscreteContinuum hybrid;

    // Can define separate directories for hybrid / discrete solutions
    hybrid.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.discreteNet.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.discreteNet.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";

    // Create build path (delete files if folder already exists)
    hybrid.discreteNet.setBuildPath(true);
    hybrid.discreteNet.loadNetwork("1Network.dat");

    // Set stack size for recursive functions (large networks)
    hybrid.discreteNet.setStackSize();

    // Use below if deadends are given in 'vesstyp' column in network file
    hybrid.discreteNet.loadDeadEnds = false;

    // Define vessel inlets / outlets (Mes. 1 example given)
    imat inOutlets = zeros<imat>(2,2);
    inOutlets(0,0) = 830;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 825;
    inOutlets(1,1) = 2;

    // Generate spatial graph for vessel classification algorithm
    hybrid.graph.generate(hybrid.discreteNet, true);

    // Run below if classifications are already assigned in 'vesstyp'
    //hybrid.graph.findTree(inOutlets);

    // Run below function if vessels haven't been classified -> run topology classification algorithm
    hybrid.graph.analyseTopology(inOutlets, hybrid.discreteNet);

    // Assign geometries to boundary nodes (based on 'vesstyp')
    hybrid.discreteNet.analyseBoundaryType();

    // Compute boundary flow and segment pressures if flow solver has not run (otherwise run below *)
    hybrid.discreteNet.bloodFlow(true);
    hybrid.discreteNet.computeSegpress();
    hybrid.discreteNet.computeBoundaryFlow();
    hybrid.discreteNet.printVisuals(false);

    // Define micro-cell build and load paths
    hybrid.cell.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";
    hybrid.cell.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    hybrid.cell.setBuildPath(false); // Create build directory / delete contents

    hybrid.graph.findLengths();

    //  Set micro-cell vessel diameter distribution(s)
    hybrid.cell.setEdgeDiamDistrib(hybrid.discreteNet.diam(find(hybrid.discreteNet.vesstyp == 2)));
    hybrid.cell.setEdgeLengthDistrib(hybrid.discreteNet.elseg(find(hybrid.discreteNet.vesstyp == 2)));
    hybrid.cell.setEucLengthDistrib(hybrid.graph.lseg(find(hybrid.graph.vesstyp == 2)));

    // Compute micro-cell hydraulic conductivity
    hybrid.cell.computeConductivity("hexCell2D");

    // Run discrete-continuum (hybrid) model
    hybrid.runHybrid();

}
