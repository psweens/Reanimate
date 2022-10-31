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
    //omp_set_num_threads(48);

    // Define vasculature
    Vasculature network;
    // Define and set load / build paths
    network.setLoadPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Load_Data/");
    network.setBuildPath("/Users/sweene01/Dropbox/Code/C++/Reanimate/Build_Data/",true);

    //network.readAmira("retina_cco_voronoicap.am", "retina_cco_voronoicap.txt");

    // Load network file
    network.graphOverride = true;
    network.loadNetwork("retina_cco_datprep_solved.txt");

    // Set stack size - used for functions w/ recursive loops (large networks)
    network.setStackSize();

    // Solve for blood flow w/ variable viscosity and phase separation
    network.bloodFlow(true, false);

    // Output 2D visuals (false -> true if Amira file needed)
    network.printVisuals(false);

}
