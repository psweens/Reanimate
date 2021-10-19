#include "Network.hpp"
#include "spatGraph.hpp"
#include "Vasculature.hpp"
#include "DiscreteContinuum.hpp"
#include "MicroCell.hpp"
#include "omp.h"

using namespace reanimate;
using namespace std;

void example_Discrete_BloodFlow() {

    // Set number of threads for any multithreading
    omp_set_num_threads(48);

    // Define vasculature
    Vasculature network;

    // Set load and built paths
    network.loadPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/";
    network.buildPath = "/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/";

    // Func. creates build path if non-existent (files deleted if folder exists)
    network.setBuildPath(true);

    // Load network file
    network.loadNetwork("1Network.dat");

    // Set stack size - used for functions w/ recursive loops (large networks)
    network.setStackSize();

    // Solve for blood flow w/ variable viscosity and phase separation
    network.bloodFlow(true, true);

    // Output 2D visuals (false -> true if Amira file needed)
    network.printVisuals(false);

}
