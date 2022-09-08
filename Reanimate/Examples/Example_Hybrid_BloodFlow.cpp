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
    //omp_set_num_threads(48);

    // Define discrete-continuum (hybrid) model
    DiscreteContinuum hybrid;

    // Define and set build / load directories
    hybrid.setBuildPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/",true);
    hybrid.setLoadPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/");

    // Load network
    hybrid.discreteNet.loadNetwork("cNetwork.txt");

    // Set stack size for recursive functions (large networks)
    hybrid.discreteNet.setStackSize();

    // Use below if deadends are given in 'vesstyp' column in network file
    hybrid.discreteNet.loadDeadEnds = false;

    // Define vessel inlets / outlets
    imat inOutlets = zeros<imat>(4,2);
    inOutlets(0,0) = 601;
    inOutlets(0,1) = 1;
    inOutlets(1,0) = 598;
    inOutlets(1,1) = 1;
    inOutlets(2,0) = 595;
    inOutlets(2,1) = 2;
    inOutlets(3,0) = 576;
    inOutlets(3,1) = 2;

    // Generate spatial graph for vessel classification algorithm
    hybrid.graph.generate(hybrid.discreteNet, true);

    // Run below if classifications are already assigned in 'vesstyp'
    hybrid.graph.findTree(inOutlets);

    // Run below function if vessels haven't been classified -> run topology classification algorithm
    //hybrid.graph.analyseTopology(inOutlets, hybrid.discreteNet);

    // Assign geometries to boundary nodes (based on 'vesstyp')
    hybrid.discreteNet.analyseBoundaryType();

    // Compute boundary flow and segment pressures if flow solver has not run (otherwise run below *)
    hybrid.discreteNet.bloodFlow(true);
    hybrid.discreteNet.computeSegpress();
    hybrid.discreteNet.computeBoundaryFlow();
    hybrid.discreteNet.printVisuals(false);

    //  Set micro-cell vessel diameter distribution(s)
    uvec indices = find_unique(hybrid.discreteNet.edgeLabels);
    ivec evesstyp = hybrid.discreteNet.vesstyp(indices);
    vec ediam = hybrid.discreteNet.ediam(indices);
    vec elseg = hybrid.discreteNet.elseg(indices);
    ediam = ediam(find(evesstyp == 2));
    elseg = elseg(find(evesstyp == 2));

    hybrid.cell.setEdgeDiamDistrib(ediam);
    hybrid.cell.setEdgeLengthDistrib(elseg);
    hybrid.graph.findLengths();
    hybrid.cell.setEucLengthDistrib(hybrid.graph.lseg(find(hybrid.graph.vesstyp == 2)));

    // Compute micro-cell hydraulic conductivity
    hybrid.cell.computeConductivity("crossCell3D", 10);

    // Run discrete-continuum (hybrid) model
    hybrid.runHybrid();

}
